//! Polytomy resolution for the timetree refinement loop.
//!
//! A polytomy is resolved by sampling a coalescent history for its children, conditioned on
//! the substitutions mapped to their branches: a lineage may only merge once every
//! substitution on its branch has been placed in time. The sweep runs backwards from the most
//! recent child toward the parent and stops when the parent's time is reached, so a polytomy
//! with too little time above it is left partly or wholly unresolved.
//!
//! [`sweep`] holds the simulation, [`apply`] the graph surgery. This module wires them to the
//! graph, collects the per-child inputs, and keeps the pre/post topology-change bookkeeping
//! the refinement loop depends on.
//!
//! Replaces the greedy pairwise merger v1 previously used (v0's `_poly`), which v0 itself
//! deprecates as unsuitable for large polytomies. Design and the v0 divergences are recorded
//! in `kb/proposals/timetree-stochastic-polytomy-resolution.md`.

pub mod apply;
pub mod sweep;

#[cfg(test)]
mod __tests__;

use crate::optimize::topology::polytomy_nodes::find_polytomy_nodes;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::timetree::optimization::polytomy::apply::{ChildRef, apply_plan};
use crate::timetree::optimization::polytomy::sweep::{Lineage, simulate_subtree};
use crate::timetree::timetree_state::TimetreeState;
use eyre::Report;
use log::debug;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::{record_merge, remove_node_if_trivial, trivial_node_branch_lengths};
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_utils::make_error;

/// Validate state required to rebuild inference after a topology change.
///
/// This runs before relaxed-clock estimation and polytomy resolution so invalid
/// input leaves the complete previous inference state untouched.
pub fn validate_tree_before_topology_change(graph: &Graph, state: &TimetreeState) -> Result<(), Report> {
  for node in graph.get_nodes() {
    if node.is_leaf() {
      continue;
    }
    let Some(time) = state.node(node.key()).time else {
      return make_error!(
        "Topology rebuild requires an inferred time for every internal node, but node {:?} has none",
        node.key()
      );
    };
    if !time.is_finite() {
      return make_error!(
        "Topology rebuild requires a finite inferred time for every internal node, but node {:?} has {time}",
        node.key()
      );
    }
  }
  Ok(())
}

/// Resolve every multifurcation in the tree by stochastic coalescent sampling.
///
/// The two rates driving the sweep are supplied by the caller, so this module holds no policy
/// about where either comes from:
///
/// - `mutation_rate` is the expected number of substitutions per unit time across the whole
///   alignment (`clock_rate * total_length`).
/// - `merger_rate` is the per-branch coalescent merger rate $\kappa(t)$ at a calendar time.
///
/// `total_length` is the summed alignment length, used only to estimate a branch's
/// substitution count when the reconstructed one cannot be read; see [`edge_mutation_count`].
///
/// Returns the number of internal nodes created. Output depends on `rng`: the same tree and
/// the same generator state produce the same topology, and different states do not.
pub fn resolve_polytomies(
  graph: &mut Graph,
  partitions: &[PartitionTimetree],
  mutation_rate: f64,
  total_length: usize,
  merger_rate: &PiecewiseConstantFn,
  rng: &mut dyn rand::RngCore,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  state: &mut TimetreeState,
) -> Result<usize, Report> {
  let polytomy_keys = find_polytomy_nodes(graph);
  if polytomy_keys.is_empty() {
    debug!("No polytomies to resolve");
    return Ok(0);
  }

  let mut total_created = 0;
  let mut topology_validated = false;

  for node_key in polytomy_keys {
    let created = resolve_single_polytomy(
      graph,
      partitions,
      node_key,
      mutation_rate,
      total_length,
      merger_rate,
      rng,
      &mut topology_validated,
      branch_lengths,
      state,
    )?;
    total_created += created;
  }

  let obsolete_count = remove_single_child_nodes(graph, branch_lengths)?;
  if obsolete_count > 0 {
    debug!("Removed {obsolete_count} obsolete single-child nodes");
  }

  graph.build()?;

  if total_created > 0 {
    debug!("Polytomy resolution introduced {total_created} new nodes");
  } else {
    debug!("Polytomies found but the sampled histories resolved none of them");
  }

  Ok(total_created)
}

#[allow(
  clippy::too_many_arguments,
  reason = "one call site; splitting would only shuffle the arguments"
)]
fn resolve_single_polytomy(
  graph: &mut Graph,
  partitions: &[PartitionTimetree],
  node_key: GraphNodeKey,
  mutation_rate: f64,
  total_length: usize,
  merger_rate: &PiecewiseConstantFn,
  rng: &mut dyn rand::RngCore,
  topology_validated: &mut bool,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  state: &mut TimetreeState,
) -> Result<usize, Report> {
  let parent_time = inferred_time(state, node_key)?;

  let children = collect_children(graph, partitions, node_key, total_length, branch_lengths, state)?;
  if children.len() < 3 {
    return Ok(0);
  }

  let lineages: Vec<Lineage> = children
    .iter()
    .map(|child| Lineage {
      time: child.time,
      mutations: child.mutations,
    })
    .collect();

  let plan = simulate_subtree(&lineages, parent_time, mutation_rate, merger_rate, rng)?;
  if plan.mergers.is_empty() {
    debug!(
      "Polytomy at node {node_key}: {} children, sampled history merged none",
      children.len()
    );
    return Ok(0);
  }

  if !*topology_validated {
    validate_tree_before_topology_change(graph, state)?;
    *topology_validated = true;
  }

  let child_refs: Vec<ChildRef> = children
    .iter()
    .map(|child| ChildRef {
      node_key: child.node_key,
      edge_key: child.edge_key,
      time: child.time,
    })
    .collect();

  let created = apply_plan(graph, node_key, parent_time, &child_refs, &plan, branch_lengths, state)?;

  debug!(
    "Polytomy at node {node_key}: {} children -> {} children, created {created} nodes",
    children.len(),
    plan.roots.len()
  );

  Ok(created)
}

/// One child of a polytomy, with everything the sweep and the surgery need.
struct ChildInfo {
  node_key: GraphNodeKey,
  edge_key: GraphEdgeKey,
  time: f64,
  mutations: u32,
}

fn collect_children(
  graph: &Graph,
  partitions: &[PartitionTimetree],
  node_key: GraphNodeKey,
  total_length: usize,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  state: &TimetreeState,
) -> Result<Vec<ChildInfo>, Report> {
  let edge_keys = {
    let node = graph.get_node(node_key).expect("Node must exist");
    node.outbound().to_vec()
  };

  edge_keys
    .into_iter()
    .map(|edge_key| {
      let edge = graph.get_edge(edge_key).expect("Edge must exist");
      let child_key = edge.target();

      let time = inferred_time(state, child_key)?;
      let mutation_length = branch_lengths.get(&edge_key).copied().flatten();

      Ok(ChildInfo {
        node_key: child_key,
        edge_key,
        time,
        mutations: edge_mutation_count(graph, partitions, edge_key, mutation_length, total_length),
      })
    })
    .collect()
}

/// Substitutions mapped to one branch.
///
/// Prefers the reconstructed substitution list, which is exact. That list is repopulated by
/// the marginal forward pass and is unavailable before the first one has run, so fall back to
/// v0's estimate from the observed branch length -- `round(mutation_length * L)` -- when it
/// cannot be read.
fn edge_mutation_count(
  graph: &Graph,
  partitions: &[PartitionTimetree],
  edge_key: GraphEdgeKey,
  mutation_length: Option<f64>,
  total_length: usize,
) -> u32 {
  let exact: Option<usize> = partitions
    .iter()
    .map(|partition| partition.edge_subs(graph, edge_key).ok().map(|subs| subs.len()))
    .sum();

  let count = exact.unwrap_or_else(|| {
    let estimate = mutation_length.unwrap_or(0.0) * total_length as f64;
    if estimate.is_finite() && estimate > 0.0 {
      estimate.round() as usize
    } else {
      0
    }
  });

  u32::try_from(count).unwrap_or(u32::MAX)
}

fn inferred_time(state: &TimetreeState, node_key: GraphNodeKey) -> Result<f64, Report> {
  let Some(time) = state.node(node_key).time else {
    return make_error!("Polytomy resolution requires an inferred time for node {node_key}, but it has none");
  };
  if !time.is_finite() {
    return make_error!("Polytomy resolution requires a finite inferred time for node {node_key}, but it has {time}");
  }
  Ok(time)
}

/// Remove single-child internal nodes left by polytomy resolution.
///
/// Uses `remove_node_if_trivial` which properly sums branch lengths into the
/// merged edge and calls `graph.build()` per removal.
fn remove_single_child_nodes(
  graph: &mut Graph,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<usize, Report> {
  let mut removed_count = 0;

  loop {
    let obsolete_key = graph.get_nodes().into_iter().find_map(|node| {
      let is_trivial = node.inbound().len() == 1 && node.outbound().len() == 1;
      is_trivial.then_some(node.key())
    });

    let Some(node_key) = obsolete_key else {
      break;
    };

    let (parent_branch, child_branch) = trivial_node_branch_lengths(graph, node_key, branch_lengths);
    if let Some(merge) = remove_node_if_trivial(graph, node_key, parent_branch, child_branch)? {
      record_merge(branch_lengths, &merge);
      removed_count += 1;
    }
  }

  Ok(removed_count)
}

/// Rebuild inference state after tree topology changes.
///
/// Internal nodes receive point constraints at their previously inferred times or
/// the times assigned during polytomy resolution. Leaf nodes retain their observed
/// date constraints and exclusion flags. Edge distributions, messages, and
/// topology-dependent rate state are reset unconditionally. Branch lengths and
/// time lengths remain valid inputs for the next inference pass.
pub fn prepare_tree_after_topology_change(graph: &Graph, state: &mut TimetreeState) -> Result<(), Report> {
  validate_tree_before_topology_change(graph, state)?;

  for node in graph.get_nodes() {
    if node.is_leaf() {
      continue;
    }
    let key = node.key();
    let Some(time) = state.node(key).time else {
      return make_error!(
        "Topology rebuild requires an inferred time for every internal node, but node {key:?} has none"
      );
    };
    // Negative-log ordinate `0` is the `NegLog` multiplicative identity (probability 1). The point
    // distribution lives in the threaded date state, the home the passes read.
    let distribution = Arc::new(Distribution::point(time, 0.0));
    state.node_mut(key).time_distribution = Some(distribution);
  }

  // Fields whose meaning depends on the previous edge topology (branch-length distribution, backward
  // message, relaxed-clock rate multiplier) live in the threaded date state and are reset there by
  // `TimetreeState::reset_date_edges_for_topology_change`. The observed branch length and inferred
  // time length are kept: both seed the next pass.

  Ok(())
}
