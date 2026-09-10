use crate::clock::clock_graph::GraphClock;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_output::write_clock_model;
use crate::clock::clock_regression::ClockParams;
use crate::clock::clock_state::ClockState;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, OptimizationMethod};
use crate::clock::pipeline::{self, ClockInput, ClockPipelineParams};
use crate::clock::rtt::{ClockRegressionResult, write_clock_regression_result_csv};
use crate::commands::clock::args::{BranchSplitArgs, TreetimeClockArgs};
use crate::commands::shared::output::OutputSelection;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::shared::tree_output::write_clock_tree_outputs;
use crate::make_error;
use crate::make_report;
use eyre::{Report, WrapErr};
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_graph::value_maps::node_names;
use treetime_io::dates_csv::read_dates;
use treetime_io::nwk::nwk_read_file;

#[derive(serde::Serialize)]
pub struct ClockGraphData {
  pub clock_model: ClockModel,
  pub regression_results: Vec<ClockRegressionResult>,
}

impl ClockGraphData {
  pub fn new(clock_model: ClockModel, regression_results: Vec<ClockRegressionResult>) -> Self {
    Self {
      clock_model,
      regression_results,
    }
  }
}

/// Per-node clock output as a value.
///
/// Holds the durable per-node results the clock output writers consume: the estimated `time`
/// (numerical date), the cumulative divergence `div`, and the two exclusion flags. `name` is carried
/// for writers that key by name. The values are gathered from the tree after clock estimation and
/// rerooting complete, so the map is keyed by the final (post-reroot) node set.
#[derive(Debug, Clone, Serialize)]
pub struct ClockNodeOut {
  pub name: Option<String>,
  pub div: f64,
  pub time: Option<f64>,
  pub is_outlier: bool,
  pub bad_branch: bool,
}

/// Per-edge output as a value: the branch length the output writers read.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct EdgeOut {
  pub branch_length: Option<f64>,
}

#[derive(serde::Serialize)]
pub struct ClockResult {
  #[serde(skip)]
  pub graph: GraphClock<ClockGraphData>,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, ClockNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
}

impl std::ops::Deref for ClockResult {
  type Target = ClockGraphData;

  fn deref(&self) -> &Self::Target {
    self.graph.data()
  }
}

/// Gather the per-node and per-edge clock outputs off the estimated tree into keyed value maps.
///
/// Runs after clock estimation, rerooting, and topology ordering, so it reads the final divergence,
/// time, and exclusion flags of the post-reroot node set. The clock passes still write these fields
/// onto the graph payloads (the reroot and best-root search read them there); this step surfaces them
/// as a standalone value the output writers consume.
fn gather_clock_outputs(
  graph: &GraphClock<ClockGraphData>,
  state: &ClockState,
) -> (BTreeMap<GraphNodeKey, ClockNodeOut>, BTreeMap<GraphEdgeKey, EdgeOut>) {
  let nodes = graph
    .get_nodes()
    .iter()
    .map(|node| {
      let node = node.read_arc();
      let key = node.key();
      // Name is an input field carried on the payload; the clock inference fields come from the
      // clock state value the pipeline routed through estimation and rerooting.
      let name = node.payload().read_arc().name.clone();
      let node_state = state.node(key);
      let out = ClockNodeOut {
        name,
        div: node_state.div,
        time: node_state.time,
        is_outlier: node_state.is_outlier,
        bad_branch: node_state.bad_branch,
      };
      (key, out)
    })
    .collect();

  let edges = graph
    .get_edges()
    .iter()
    .map(|edge| {
      let edge = edge.read_arc();
      let branch_length = edge.payload().read_arc().branch_length;
      (edge.key(), EdgeOut { branch_length })
    })
    .collect();

  (nodes, edges)
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
    clock_args.metadata(),
    &clock_args.metadata_id.metadata_delimiters,
    &clock_args.metadata_id.metadata_id_columns,
    &None,
    &clock_args.date_column.date_column,
  )
  .wrap_err("When reading dates")?;

  let resolved = clock_args.resolve_outputs()?;

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
    state,
    clock_model,
    regression_results,
  } = output;
  let mut graph = graph.map_data(ClockGraphData::new(clock_model, regression_results));
  let topology_order = clock_args
    .topology_order
    .resolve_topology_order(&graph, Some(input_order))?;
  topology_order.apply(&mut graph)?;
  progress.report("Writing output", 0.8, "");

  let (nodes, edges) = gather_clock_outputs(&graph, &state);
  let branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>> =
    edges.iter().map(|(key, edge)| (*key, edge.branch_length)).collect();

  if !resolved.tree_outputs.is_empty() {
    write_clock_tree_outputs(
      &graph,
      &nodes,
      &branch_lengths,
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
  Ok(ClockResult { graph, nodes, edges })
}

fn leaf_order<N, E, D>(graph: &Graph<N, E, D>) -> Result<Vec<String>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  let names = node_names(graph);
  graph
    .get_leaves()
    .into_iter()
    .map(|leaf| {
      let key = leaf.read_arc().key();
      names[&key]
        .clone()
        .ok_or_else(|| make_report!("Leaf node {key} has no name"))
    })
    .collect()
}
