use crate::clock::clock_graph::GraphClock;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_output::write_clock_model;
use crate::clock::clock_regression::ClockParams;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, OptimizationMethod};
use crate::clock::pipeline::{self, ClockInput, ClockPipelineParams};
use crate::clock::rtt::{ClockRegressionResult, write_clock_regression_result_csv};
use crate::commands::clock::args::{BranchSplitArgs, TreetimeClockArgs};
use crate::commands::shared::output::{CommandKind, OutputSelection};
use crate::commands::shared::tree_output::TreeOutputAdapter;
use crate::make_error;
use crate::make_report;
use eyre::{Report, WrapErr};
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_io::dates_csv::read_dates;
use treetime_io::graph::write_tree_outputs;
use treetime_io::nwk::nwk_read_file;

#[allow(clippy::manual_non_exhaustive, clippy::partial_pub_fields)] // The private unit field preserves Graph JSON's `data: null` shape.
#[derive(serde::Serialize)]
#[serde(transparent)]
pub struct ClockGraphData {
  marker: (),
  #[serde(skip)]
  pub clock_model: ClockModel,
  #[serde(skip)]
  pub regression_results: Vec<ClockRegressionResult>,
}

impl ClockGraphData {
  pub fn new(clock_model: ClockModel, regression_results: Vec<ClockRegressionResult>) -> Self {
    Self {
      marker: (),
      clock_model,
      regression_results,
    }
  }
}

#[derive(serde::Serialize)]
pub struct ClockResult {
  #[serde(skip)]
  pub graph: GraphClock<ClockGraphData>,
}

impl std::ops::Deref for ClockResult {
  type Target = ClockGraphData;

  fn deref(&self) -> &Self::Target {
    self.graph.data()
  }
}

fn branch_split_to_params(args: &BranchSplitArgs) -> BranchPointOptimizationParams {
  match args.method {
    OptimizationMethod::Grid => BranchPointOptimizationParams::grid_with(args.grid_params.clone()),
    OptimizationMethod::Brent => BranchPointOptimizationParams::brent_with(args.brent_params.clone()),
    OptimizationMethod::GoldenSection => BranchPointOptimizationParams::golden_section_with(args.golden_params.clone()),
  }
}

pub fn run_clock(
  clock_args: &TreetimeClockArgs,
  progress: &dyn crate::progress::ProgressSink,
) -> Result<ClockResult, Report> {
  progress.check_cancelled()?;
  progress.report("Reading input", 0.0, "");

  let graph: GraphClock = if let Some(tree) = &clock_args.tree {
    nwk_read_file(tree)
  } else {
    return make_error!("Tree inference is not implemented. Provide a tree file with --tree");
  }?;
  let input_order = leaf_order(&graph)?;

  let dates = read_dates(
    &clock_args.metadata,
    &clock_args.metadata_id.metadata_id_columns,
    &None,
    &clock_args.date_column.date_column,
  )
  .wrap_err("When reading dates")?;

  let selection: Vec<OutputSelection> = clock_args
    .output_selection
    .iter()
    .copied()
    .map(OutputSelection::from)
    .collect();
  let resolved = clock_args.output.resolve(
    CommandKind::Clock,
    &selection,
    &[
      (OutputSelection::ClockModel, clock_args.output_clock_model.as_deref()),
      (OutputSelection::ClockCsv, clock_args.output_clock_csv.as_deref()),
    ],
  )?;

  let clock_params = if clock_args.covariation {
    let seq_len = clock_args
      .sequence_length
      .ok_or_else(|| make_report!("--sequence-length is required when --covariation is enabled"))?
      as f64;
    let tip_slack = clock_args.tip_slack.unwrap_or(3.0);
    let overdispersion = 2.0;
    ClockParams {
      variance_factor: overdispersion / seq_len,
      variance_offset: 0.0,
      variance_offset_leaf: tip_slack * tip_slack / seq_len / seq_len,
    }
  } else {
    clock_args.clock_regression.clock_params.clone()
  };

  let params = ClockPipelineParams {
    clock_params,
    clock_filter: clock_args.clock_filter,
    keep_root: clock_args.keep_root,
    allow_negative_rate: clock_args.allow_negative_rate,
    branch_params: branch_split_to_params(&clock_args.branch_split),
    reroot_spec: clock_args.reroot.spec(),
  };

  let input = ClockInput { graph, dates };

  let output = pipeline::run(&params, input, progress)?;
  let pipeline::ClockOutput {
    graph,
    clock_model,
    regression_results,
  } = output;
  let mut graph = graph.map_data(ClockGraphData::new(clock_model, regression_results));
  let topology_order = clock_args
    .topology_order
    .resolve_topology_order(&graph, Some(input_order))?;
  topology_order.apply(&mut graph)?;
  resolved.prepare()?;

  progress.report("Writing output", 0.8, "");

  if !resolved.tree_outputs.is_empty() {
    write_tree_outputs::<TreeOutputAdapter, _, _, _>(
      &graph,
      &resolved.tree_outputs,
      &treetime_io::nwk::CommentProviders::new(),
    )?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ClockModel) {
    write_clock_model(&graph.data().clock_model, path)?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ClockCsv) {
    write_clock_regression_result_csv(&graph.data().regression_results, path, b',')?;
  }

  progress.report("Done", 1.0, "");
  Ok(ClockResult { graph })
}

fn leaf_order<N, E, D>(graph: &Graph<N, E, D>) -> Result<Vec<String>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  graph
    .get_leaves()
    .into_iter()
    .map(|leaf| {
      let leaf = leaf.read_arc();
      leaf
        .payload()
        .read_arc()
        .name()
        .map(|name| name.as_ref().to_owned())
        .ok_or_else(|| make_report!("Leaf node {} has no name", leaf.key()))
    })
    .collect()
}
