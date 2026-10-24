use crate::clock::clock_state::{ClockEdgeInput, ClockInputs, ClockNodeInput};
use crate::progress::LogSink;
use crate::progress_info;
use crate::reroot::find_root::find_root;
use crate::reroot::orchestrate::apply_root_at_edge;
use crate::reroot::params::RerootParams;
use crate::reroot::placement::root_moves;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::reroot::{RerootResult, StemRemovalInfo, remove_stem_root};

pub(crate) fn reroot_at_best_root(
  graph: &mut Graph,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  inputs: &mut ClockInputs,
  outliers: &BTreeSet<GraphNodeKey>,
  reroot_params: &RerootParams,
  names: &mut NodeNames,
  log: &dyn LogSink,
) -> Result<RerootResult, Report> {
  let split_edge = reroot_params.topology.split_edge;
  let placement = reroot_params.spec.placement();
  let mut point = find_root(graph, branch_lengths, inputs, outliers, reroot_params, log)?;
  let mut stem_removal = None;
  if root_moves(graph, branch_lengths, point, placement, split_edge)? {
    stem_removal = remove_undated_stem(graph, branch_lengths, inputs)?;
    if stem_removal.is_some() {
      point = find_root(graph, branch_lengths, inputs, outliers, reroot_params, log)?;
    }
  }
  let result = apply_root_at_edge(
    graph,
    point,
    placement,
    reroot_params.topology,
    branch_lengths,
    stem_removal,
    |key| inputs.likely_time(key).is_some(),
    names,
  )?;
  update_clock_inputs(inputs, &result);
  progress_info!(log, "Rerooted to {}", &names[result.new_root_key]);
  Ok(result)
}

pub(crate) fn update_clock_inputs(inputs: &mut ClockInputs, result: &RerootResult) {
  if let Some(stem) = &result.stem_removal {
    inputs.nodes.remove(&stem.removed_node_key);
    inputs.edges.remove(&stem.removed_edge_key);
  }
  if let Some(split) = &result.edge_split {
    inputs.edges.remove(&split.old_edge_key);
    inputs
      .edges
      .insert(split.parent_side_edge_key, ClockEdgeInput::default());
    inputs
      .edges
      .insert(split.child_side_edge_key, ClockEdgeInput::default());
    inputs.nodes.insert(split.new_node_key, ClockNodeInput::default());
  }
  if let Some(merge) = &result.edge_merge {
    inputs.nodes.remove(&merge.removed_node_key);
    inputs.edges.remove(&merge.parent_edge_key);
    inputs.edges.remove(&merge.child_edge_key);
    inputs.edges.insert(merge.merged_edge_key, ClockEdgeInput::default());
  }
}

fn remove_undated_stem(
  graph: &mut Graph,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  inputs: &mut ClockInputs,
) -> Result<Option<StemRemovalInfo>, Report> {
  let root_key = graph.root_key()?;
  if inputs.likely_time(root_key).is_some() {
    return Ok(None);
  }
  let Some(stem) = remove_stem_root(graph, root_key)? else {
    return Ok(None);
  };
  branch_lengths.remove(&stem.removed_edge_key);
  inputs.nodes.remove(&stem.removed_node_key);
  inputs.edges.remove(&stem.removed_edge_key);
  Ok(Some(stem))
}
