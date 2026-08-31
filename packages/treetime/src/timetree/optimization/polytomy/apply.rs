//! Realise a [`SubtreePlan`] on the graph.
//!
//! Kept apart from [`super::sweep`] so the simulation stays free of graph state. Every
//! mutation of the tree for one polytomy happens here, in a single pass over the plan.

use crate::partition::timetree::partition::GraphTimetree;
use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
use crate::timetree::optimization::polytomy::sweep::SubtreePlan;
use eyre::Report;
use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};
use treetime_graph::node::GraphNodeKey;
use treetime_utils::{make_internal_error, make_internal_report};

/// One child of the polytomy, paired with the graph objects the plan refers to by index.
#[derive(Clone, Copy, Debug)]
pub struct ChildRef {
  pub node_key: GraphNodeKey,
  pub edge_key: GraphEdgeKey,
  pub time: f64,
}

/// Apply every merger in `plan`, then reattach the surviving lineages to `parent_key`.
///
/// Returns the number of nodes created, which equals `plan.mergers.len()`.
///
/// Mergers are processed in order. Each creates a node and pulls its two lineages beneath it;
/// because a merger may only reference earlier mergers, the node a lineage needs always
/// exists by the time it is referenced.
///
/// Original children are relocated with [`treetime_graph::graph::Graph::reparent_edge`], which
/// preserves the edge key and payload. That matters here: the child's `branch_length` is the
/// observed mutation length, which `prepare_tree_after_topology_change` deliberately keeps to
/// seed the next inference pass, and partition state is keyed by edge key, so a fresh key
/// would have its entry rebuilt as empty. Only `time_length` is rewritten, since only the
/// parent moved.
pub fn apply_plan(
  graph: &mut GraphTimetree,
  parent_key: GraphNodeKey,
  parent_time: f64,
  children: &[ChildRef],
  plan: &SubtreePlan,
) -> Result<usize, Report> {
  // Validate the complete forest before graph mutation so plan application is atomic.
  let times = validate_plan(parent_time, children, plan)?;

  let mut merger_nodes: Vec<GraphNodeKey> = Vec::with_capacity(plan.mergers.len());

  for merger in &plan.mergers {
    let new_node_key = graph.add_node(NodeTimetree {
      time: Some(merger.time),
      ..NodeTimetree::default()
    });

    for lineage in [merger.left, merger.right] {
      attach(
        graph,
        children,
        &merger_nodes,
        &times,
        lineage,
        new_node_key,
        merger.time,
      )?;
    }

    merger_nodes.push(new_node_key);
  }

  for &lineage in &plan.roots {
    attach(graph, children, &merger_nodes, &times, lineage, parent_key, parent_time)?;
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

/// Place one lineage under `new_parent_key`, setting the connecting edge's `time_length`.
fn attach(
  graph: &mut GraphTimetree,
  children: &[ChildRef],
  merger_nodes: &[GraphNodeKey],
  times: &[f64],
  lineage: usize,
  new_parent_key: GraphNodeKey,
  new_parent_time: f64,
) -> Result<(), Report> {
  let Some(&lineage_time) = times.get(lineage) else {
    return make_internal_error!("Polytomy plan referenced unknown lineage {lineage}");
  };
  let time_length = lineage_time - new_parent_time;

  if let Some(child) = children.get(lineage) {
    // An original child: relocate its existing edge, keeping key and payload.
    graph.reparent_edge(child.edge_key, new_parent_key)?;
    let edge = graph
      .get_edge(child.edge_key)
      .ok_or_else(|| make_internal_report!("Edge {} vanished while applying polytomy plan", child.edge_key))?;
    edge.write_arc().payload().write_arc().time_length = Some(time_length);
  } else {
    // A node the sweep created: it has no parent edge yet.
    let Some(&node_key) = merger_nodes.get(lineage - children.len()) else {
      return make_internal_error!(
        "Polytomy plan referenced merger node {lineage} before it was created; mergers must only reference earlier mergers"
      );
    };
    let mut payload = EdgeTimetree {
      time_length: Some(time_length),
      ..EdgeTimetree::default()
    };
    // The sweep only merges lineages that have placed every substitution, so the branch
    // above a merger node carries none.
    payload.set_branch_length(Some(0.0));
    graph.add_edge(new_parent_key, node_key, payload)?;
  }

  Ok(())
}
