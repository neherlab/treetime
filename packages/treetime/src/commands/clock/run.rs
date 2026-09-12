use crate::clock::clock_model::ClockModel;
use crate::clock::clock_output::write_clock_model;
use crate::clock::clock_regression::ClockParams;
use crate::clock::clock_state::{ClockInputs, ClockState};
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
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::dates_csv::read_dates;
use treetime_io::nwk::{NwkParse, nwk_read_file};

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
  pub graph: Graph,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, ClockNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
  #[serde(skip)]
  pub clock_model: ClockModel,
  #[serde(skip)]
  pub regression_results: Vec<ClockRegressionResult>,
}

/// Gather the per-node and per-edge clock outputs into keyed value maps the output writers consume.
///
/// Runs after clock estimation, rerooting, and topology ordering. Each node's divergence, time, and
/// exclusion flags come from the `ClockState` value the pipeline routed through estimation and
/// rerooting; each node's name from the post-reroot `names` map; each edge's branch length from the
/// `branch_lengths` map. The maps are keyed by the final
/// (post-reroot) node and edge set.
fn gather_clock_outputs(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &ClockState,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> (BTreeMap<GraphNodeKey, ClockNodeOut>, BTreeMap<GraphEdgeKey, EdgeOut>) {
  let nodes = graph
    .get_nodes()
    .iter()
    .map(|node| {
      let key = node.read_arc().key();
      // Name comes from the post-reroot name map the pipeline returns; the fitted clock results
      // (divergence, outlier flag) come from the clock state value, while the observed date and
      // bad-branch flag come from the clock inputs value the pipeline routed through estimation and
      // rerooting.
      let name = names[&key].clone();
      let node_state = state.node(key);
      let node_input = inputs.node(key);
      let out = ClockNodeOut {
        name,
        div: node_state.div,
        time: node_input.time,
        is_outlier: node_state.is_outlier,
        bad_branch: node_input.bad_branch,
      };
      (key, out)
    })
    .collect();

  let edges = graph
    .get_edges()
    .iter()
    .map(|edge| {
      let key = edge.read_arc().key();
      (
        key,
        EdgeOut {
          branch_length: branch_lengths[&key],
        },
      )
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

  let NwkParse {
    graph,
    names,
    branch_lengths,
    ..
  } = if let Some(tree) = &clock_args.tree {
    nwk_read_file(tree)
  } else {
    return make_error!("Tree inference is not implemented. Provide a tree file with --tree");
  }?;
  let input_order = leaf_order(&graph, &names)?;

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

  let input = ClockInput {
    graph,
    dates,
    branch_lengths,
  };

  let output = pipeline::run(&params, input, &names, progress)?;
  let pipeline::ClockOutput {
    mut graph,
    inputs,
    state,
    clock_model,
    regression_results,
    names,
    branch_lengths,
  } = output;
  let topology_order = clock_args
    .topology_order
    .resolve_topology_order(&graph, &names, Some(input_order))?;
  topology_order.apply(&mut graph, &names, &branch_lengths)?;
  progress.report("Writing output", 0.8, "");

  // The pipeline's post-reroot name and branch-length maps carry the final tree's values; topology
  // ordering only permutes children, so the maps still match after `apply`. Both drive the gather
  // and the tree-output writers below.
  let (nodes, edges) = gather_clock_outputs(&graph, &inputs, &state, &names, &branch_lengths);

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
    write_clock_model(&clock_model, path)?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ClockCsv) {
    write_clock_regression_result_csv(&regression_results, path, b',')?;
  }

  progress.report("Done", 1.0, "");
  Ok(ClockResult {
    graph,
    nodes,
    edges,
    clock_model,
    regression_results,
  })
}

fn leaf_order(graph: &Graph, names: &BTreeMap<GraphNodeKey, Option<String>>) -> Result<Vec<String>, Report> {
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
