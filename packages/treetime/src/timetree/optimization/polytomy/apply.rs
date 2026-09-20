use crate::timetree::optimization::polytomy::sweep::SubtreePlan;
use crate::timetree::timetree_state::{DateNodeState, TimetreeState};
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::make_internal_error;

#[derive(Clone, Copy, Debug)]
pub struct ChildRef {
  pub node_key: GraphNodeKey,
  pub edge_key: GraphEdgeKey,
  pub time: f64,
}

pub fn apply_plan(
  graph: &mut Graph,
  parent_key: GraphNodeKey,
  parent_time: f64,
  children: &[ChildRef],
  plan: &SubtreePlan,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  state: &mut TimetreeState,
) -> Result<usize, Report> {
  let times = validate_plan(parent_time, children, plan)?;

  let mut merger_nodes: Vec<GraphNodeKey> = Vec::with_capacity(plan.mergers.len());

  for merger in &plan.mergers {
    let new_node_key = graph.add_node();
    state.nodes.insert(
      new_node_key,
      DateNodeState {
        time: Some(merger.time),
        ..DateNodeState::default()
      },
    );

    for lineage in [merger.left, merger.right] {
      attach(
        graph,
        children,
        &merger_nodes,
        &times,
        lineage,
        new_node_key,
        merger.time,
        branch_lengths,
        state,
      )?;
    }

    merger_nodes.push(new_node_key);
  }

  for &lineage in &plan.roots {
    attach(
      graph,
      children,
      &merger_nodes,
      &times,
      lineage,
      parent_key,
      parent_time,
      branch_lengths,
      state,
    )?;
  }

  Ok(plan.mergers.len())
}

fn validate_plan(parent_time: f64, children: &[ChildRef], plan: &SubtreePlan) -> Result<Vec<f64>, Report> {
  if !parent_time.is_finite() {
    return make_internal_error!("Polytomy plan parent time must be finite, got {parent_time}");
  }

  let n_children = children.len();
  let n_lineages = n_children + plan.mergers.len();
  let mut times = Vec::with_capacity(n_lineages);
  for (index, child) in children.iter().enumerate() {
    if !child.time.is_finite() || child.time < parent_time {
      return make_internal_error!(
        "Polytomy plan child {index} time must be finite and not older than its parent at {parent_time:.6e}, got {:.6e}",
        child.time
      );
    }
    times.push(child.time);
  }

  let mut consumed = vec![false; n_lineages];
  for (index, merger) in plan.mergers.iter().enumerate() {
    if !merger.time.is_finite() || merger.time <= parent_time {
      return make_internal_error!(
        "Polytomy plan merger {index} time must be finite and more recent than its parent at {parent_time:.6e}, got {:.6e}",
        merger.time
      );
    }

    let merger_id = n_children + index;
    for lineage in [merger.left, merger.right] {
      if lineage >= merger_id {
        return make_internal_error!("Polytomy plan merger {index} referenced lineage {lineage} before it was created");
      }
      if consumed[lineage] {
        return make_internal_error!("Polytomy plan consumed lineage {lineage} more than once");
      }
      if times[lineage] < merger.time {
        return make_internal_error!(
          "Polytomy plan merger {index} at {:.6e} is more recent than lineage {lineage} at {:.6e}",
          merger.time,
          times[lineage]
        );
      }
      consumed[lineage] = true;
    }
    times.push(merger.time);
  }

  for &lineage in &plan.roots {
    if lineage >= n_lineages {
      return make_internal_error!("Polytomy plan referenced unknown root lineage {lineage}");
    }
    if consumed[lineage] {
      return make_internal_error!("Polytomy plan consumed lineage {lineage} more than once");
    }
    if times[lineage] < parent_time {
      return make_internal_error!(
        "Polytomy plan root lineage {lineage} at {:.6e} is older than its parent at {parent_time:.6e}",
        times[lineage]
      );
    }
    consumed[lineage] = true;
  }

  if let Some(lineage) = consumed.iter().position(|consumed| !consumed) {
    return make_internal_error!("Polytomy plan omitted lineage {lineage}");
  }

  Ok(times)
}

fn attach(
  graph: &mut Graph,
  children: &[ChildRef],
  merger_nodes: &[GraphNodeKey],
  times: &[f64],
  lineage: usize,
  new_parent_key: GraphNodeKey,
  new_parent_time: f64,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  state: &mut TimetreeState,
) -> Result<(), Report> {
  let Some(&lineage_time) = times.get(lineage) else {
    return make_internal_error!("Polytomy plan referenced unknown lineage {lineage}");
  };
  let time_length = lineage_time - new_parent_time;

  if let Some(child) = children.get(lineage) {
    graph.reparent_edge(child.edge_key, new_parent_key)?;
    state.edges.entry(child.edge_key).or_default().time_length = Some(time_length);
  } else {
    let Some(&node_key) = merger_nodes.get(lineage - children.len()) else {
      return make_internal_error!(
        "Polytomy plan referenced merger node {lineage} before it was created; mergers must only reference earlier mergers"
      );
    };
    let new_edge_key = graph.add_edge(new_parent_key, node_key)?;
    branch_lengths.insert(new_edge_key, Some(0.0));
    state.edges.entry(new_edge_key).or_default().time_length = Some(time_length);
  }

  Ok(())
}
