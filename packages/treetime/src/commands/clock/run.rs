use crate::clock::clock_graph::GraphClock;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_output::write_clock_model;
use crate::clock::clock_regression::ClockParams;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, OptimizationMethod};
use crate::clock::pipeline::{self, ClockInput, ClockPipelineParams};
use crate::clock::rtt::{ClockRegressionResult, write_clock_regression_result_csv};
use crate::commands::clock::args::{BranchSplitArgs, TreetimeClockArgs};
use crate::make_error;
use crate::make_report;
use eyre::{Report, WrapErr};
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_io::dates_csv::read_dates;
use treetime_io::graph::write_graph_files_with_options;
use treetime_io::nwk::nwk_read_file;

#[derive(Debug, serde::Serialize)]
pub struct ClockResult {
  #[serde(skip)]
  pub graph: GraphClock,
  pub clock_model: ClockModel,
  pub regression_results: Vec<ClockRegressionResult>,
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

  progress.report("Writing output", 0.8, "");
  let outdir = &clock_args.output.outdir;
  let graph_options = clock_args
    .output
    .graph_write_options_with_input_order(&output.graph, input_order)?;
  write_graph_files_with_options(
    outdir,
    "rerooted",
    &output.graph,
    &treetime_io::nwk::CommentProviders::new(),
    &graph_options,
  )?;
  write_clock_model(&output.clock_model, &outdir.join("clock_model"))?;
  write_clock_regression_result_csv(&output.regression_results, outdir.join("clock.csv"), b',')?;

  progress.report("Done", 1.0, "");
  Ok(ClockResult {
    graph: output.graph,
    clock_model: output.clock_model,
    regression_results: output.regression_results,
  })
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
