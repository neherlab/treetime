use crate::optimize::topology::collapse::collapse_edge;
use crate::optimize::topology::hoist_reversions::{count_child_reversions, hoist_reverting_child};
use crate::optimize::topology::merge_shared_mutations::merge_single_polytomy;
use crate::optimize::topology::polytomy_nodes::find_polytomy_nodes;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::payload::ancestral::GraphAncestral;
use eyre::Report;
use log::debug;
use parking_lot::RwLock;
use std::collections::BTreeSet;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;

/// Resolve reversion-driven and shared-mutation polytomies across the whole tree.
///
/// Runs [`resolve_one`] over every polytomy and repeats until a full pass changes nothing.
/// Each applied move (a shared-mutation merge, a reverting-child hoist, or a helper-node
/// retirement) strictly decreases the lexicographic potential (total fitch mutation count,
/// then node count), so the fixpoint is reached in a bounded number of rounds. Retirement can
/// turn a former helper's parent into a new polytomy, which the outer loop then picks up.
///
/// Sparse-only: dense partitions carry no per-edge mutation lists, so the routine is inert
/// when no sparse partition is present, matching [`merge_shared_mutation_branches`].
///
/// Returns the number of polytomies whose local structure changed, summed across rounds; a
/// non-zero result means the caller must rebuild the graph and reassign node names.
///
/// [`merge_shared_mutation_branches`]: crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches
pub fn resolve_polytomies(
  graph: &mut GraphAncestral,
  sparse: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense: &[Arc<RwLock<PartitionMarginalDense>>],
) -> Result<usize, Report> {
  if sparse.is_empty() {
    return Ok(0);
  }

  let mut total_changed = 0;
  loop {
    let polytomy_keys = find_polytomy_nodes(graph);
    let mut round_changed = 0;
    for node_key in polytomy_keys {
      if resolve_one(graph, sparse, dense, node_key)? {
        round_changed += 1;
      }
    }
    if round_changed == 0 {
      break;
    }
    total_changed += round_changed;
  }

  if total_changed > 0 {
    debug!("Resolved {total_changed} polytomies via merge/hoist/retire");
  }

  Ok(total_changed)
}

/// Apply the merge -> hoist -> retire routine at one polytomy until it stops changing.
///
/// 1. Merge siblings that share substitutions under helper nodes
///    ([`merge_single_polytomy`]). This canonicalizes several children reverting the same
///    position into a single reverting child.
/// 2. Hoist the child with the largest reversion count against the node's parent edge
///    ([`hoist_reverting_child`]), removing one mutation per reverted position. One child per
///    round keeps the move greedy and deterministic (ties broken by edge key).
/// 3. Retire helper nodes: collapse mutation-free edges whose target was created during this
///    invocation ([`retire_created_helpers`]). This dissolves the empty helper the hoist
///    leaves behind, turning its children into genuine siblings.
///
/// The routine skips the root (no parent edge to revert). Nodes created before this call are
/// recorded so retirement never dissolves a pre-existing subtree root on a mutation-free edge.
///
/// Returns whether anything changed.
fn resolve_one(
  graph: &mut GraphAncestral,
  sparse: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense: &[Arc<RwLock<PartitionMarginalDense>>],
  node_key: GraphNodeKey,
) -> Result<bool, Report> {
  let preexisting: BTreeSet<GraphNodeKey> = graph.get_nodes().iter().map(|node| node.read_arc().key()).collect();

  let mut any_changed = false;
  loop {
    let merged = merge_single_polytomy(graph, sparse, node_key)? > 0;
    let hoisted = try_hoist_reverting_child(graph, sparse, dense, node_key)?;
    let retired = retire_created_helpers(graph, sparse, dense, &preexisting)?;

    if !(merged || hoisted || retired) {
      break;
    }
    any_changed = true;
  }

  Ok(any_changed)
}

/// Hoist the best reverting child of a node, if the node has a parent edge and any child
/// reverts one of its substitutions. Returns whether a hoist was applied.
fn try_hoist_reverting_child(
  graph: &mut GraphAncestral,
  sparse: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense: &[Arc<RwLock<PartitionMarginalDense>>],
  node_key: GraphNodeKey,
) -> Result<bool, Report> {
  let Some(parent_edge_key) = single_inbound_edge(graph, node_key) else {
    return Ok(false);
  };
  let Some(child_edge_key) = best_reverting_child(graph, sparse, node_key, parent_edge_key) else {
    return Ok(false);
  };
  hoist_reverting_child(graph, sparse, dense, parent_edge_key, child_edge_key)?;
  Ok(true)
}

/// Pick the child edge with the most reversions against the parent edge.
///
/// Chooses the largest reversion count, breaking ties by smallest edge key for determinism.
/// Returns `None` when no child reverts any parent substitution.
fn best_reverting_child(
  graph: &GraphAncestral,
  sparse: &[Arc<RwLock<PartitionMarginalSparse>>],
  node_key: GraphNodeKey,
  parent_edge_key: GraphEdgeKey,
) -> Option<GraphEdgeKey> {
  let child_edges = graph.get_node(node_key)?.read_arc().outbound().to_vec();

  let mut best: Option<(usize, GraphEdgeKey)> = None;
  for child_edge_key in child_edges {
    let reversions = count_child_reversions(sparse, parent_edge_key, child_edge_key);
    if reversions == 0 {
      continue;
    }
    let better = match best {
      Some((best_reversions, best_key)) => {
        best_reversions > reversions || (best_reversions == reversions && best_key <= child_edge_key)
      },
      None => false,
    };
    if !better {
      best = Some((reversions, child_edge_key));
    }
  }

  best.map(|(_, key)| key)
}

/// Collapse mutation-free edges whose target was created during the current [`resolve_one`].
///
/// The hoist leaves an empty edge to the reverting child exactly when that child consisted of
/// nothing but reversions, which is the normal case after a merge. Collapsing it retires the
/// helper node and makes its children genuine siblings.
///
/// Restricting the target to nodes created in this invocation is essential: a merge's helper
/// edges to pre-existing subtree roots are frequently mutation-free too, and collapsing those
/// would flatten topology the input asserted and discard its branch lengths. Within this scope
/// the branch-length-zero requirement of the loop's zero-optimal collapse is relaxed, because
/// these edges were synthesized moments earlier and carry no optimizer decision to override.
///
/// Returns whether any edge was retired.
fn retire_created_helpers(
  graph: &mut GraphAncestral,
  sparse: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense: &[Arc<RwLock<PartitionMarginalDense>>],
  preexisting: &BTreeSet<GraphNodeKey>,
) -> Result<bool, Report> {
  let mut retired = false;
  loop {
    let candidate = graph.get_edges().iter().find_map(|edge| {
      let edge = edge.read_arc();
      let target_key = edge.target();
      if preexisting.contains(&target_key) {
        return None;
      }
      let target_is_leaf = graph.get_node(target_key).is_some_and(|node| node.read_arc().is_leaf());
      if target_is_leaf {
        return None;
      }
      let edge_key = edge.key();
      let mutation_free = sparse.iter().all(|partition| {
        let partition = partition.read_arc();
        match partition.edges.get(&edge_key) {
          Some(edge_data) => edge_data.fitch_subs().is_empty() && edge_data.indels.is_empty(),
          None => true,
        }
      });
      mutation_free.then_some(edge_key)
    });

    match candidate {
      Some(edge_key) => {
        collapse_edge(graph, sparse, dense, edge_key)?;
        retired = true;
      },
      None => break,
    }
  }
  Ok(retired)
}

/// The single parent edge of a node, or `None` for the root.
fn single_inbound_edge(graph: &GraphAncestral, node_key: GraphNodeKey) -> Option<GraphEdgeKey> {
  let node = graph.get_node(node_key)?;
  let node = node.read_arc();
  match node.inbound() {
    [edge_key] => Some(*edge_key),
    _ => None,
  }
}
