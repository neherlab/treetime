use crate::hacks::fix_branch_length::fix_branch_length;
use crate::partition::indexed_pass::{IndexedPassDependencies, IndexedPassSlot};
use crate::partition::marginal::shared::normalize::{forward_log_lh_add_normalization, forward_log_lh_remove_child};
use crate::partition::marginal::sparse::message::{
  combine_messages, normalize_1d_inplace, propagate_raw, propagate_raw_per_site,
};
use crate::partition::marginal::sparse::reconstruct::reconstruct_map_seq;
use crate::partition::storage::sparse::{SparseEdgePartition, SparseNodePartition, SparseSeqDistribution, VarPos};
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use std::collections::BTreeSet;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_primitives::{AsciiChar, LogLh};
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::interval::range::range_contains;

#[allow(clippy::too_many_arguments)]
pub(crate) fn process_node_forward_indexed<N, E>(
  graph: &Graph<N, E, ()>,
  alphabet: &crate::alphabet::alphabet::Alphabet,
  gtr: &crate::gtr::gtr::GTR,
  length: usize,
  root_sequence: &treetime_primitives::Seq,
  dependencies: &IndexedPassDependencies<SparseNodePartition, SparseEdgePartition>,
  slot: &mut IndexedPassSlot<SparseNodePartition, SparseEdgePartition>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  if let Some(parent_key) = slot.parent_key {
    let parent = dependencies.node(parent_key);
    let (edge_key, edge_data) = slot
      .parent_edge
      .as_mut()
      .expect("Non-root node must own its parent edge");

    edge_data.msg_to_child = compute_msg_to_child(&slot.node, parent, edge_data)?;

    let mut variable_pos = btreemap! {};
    let mut parent_state = btreemap! {};
    let mut child_state = btreemap! {};
    for mutation in edge_data.fitch_subs() {
      let current_state = mutation.qry();
      variable_pos.insert(mutation.pos(), current_state);
      parent_state.entry(mutation.pos()).or_insert_with(|| mutation.reff());
      child_state.entry(mutation.pos()).or_insert(current_state);
    }
    for (pos, profile) in &edge_data.msg_to_child.variable {
      if !range_contains(&slot.node.seq.non_char, *pos) {
        variable_pos.entry(*pos).or_insert(profile.state);
        parent_state.entry(*pos).or_insert(profile.state);
      }
    }
    for (pos, profile) in &edge_data.msg_to_parent.variable {
      variable_pos.entry(*pos).or_insert(profile.state);
      child_state.entry(*pos).or_insert(profile.state);
    }

    let branch_length = graph
      .get_edge(*edge_key)
      .expect("Indexed edge must exist in graph")
      .read_arc()
      .payload()
      .read_arc()
      .profile_branch_length()
      .unwrap_or(0.0);
    let branch_length = fix_branch_length(length, branch_length);
    let msg_from_parent = if gtr.has_site_rates() {
      propagate_raw_per_site(
        gtr,
        branch_length,
        false,
        &edge_data.msg_to_child,
        edge_data.transmission.as_deref(),
      )
    } else {
      propagate_raw(
        &gtr.expQt(branch_length),
        &edge_data.msg_to_child,
        edge_data.transmission.as_deref(),
      )
    };
    // Persist the down-message only for tips: reconstruct_node_sequence imputes missing tip states
    // from it and has no branch length to recompute the propagation. Internal nodes never need it.
    if graph.is_leaf(slot.key) {
      edge_data.msg_from_parent = msg_from_parent.clone();
    }
    let profile = combine_messages(
      &slot.node.seq.composition,
      &[msg_from_parent, edge_data.msg_to_parent.clone()],
      &variable_pos,
      &[parent_state, child_state],
      alphabet,
      None,
    )?;
    slot.node.profile = profile;

    // Internal nodes derive their MAP sequence from the parent. Leaves keep their observed input:
    // deriving a leaf from the parent discards observed states the leaf shares with the parent under
    // Fitch compression, which corrupts both the emitted tip sequence and the leaf-edge `ml_subs`
    // (dropping real parent->tip substitutions). The dense backend likewise preserves the observed
    // leaf sequence through its forward pass.
    if !graph.is_leaf(slot.key) && !parent.seq.sequence.is_empty() {
      slot.node.seq.sequence = reconstruct_map_seq(&parent.seq.sequence, Some(edge_data), &slot.node, alphabet);
    }
    edge_data.set_ml_subs(compute_ml_subs_for_nodes(alphabet, parent, &slot.node, edge_data)?);
  } else if slot.node.seq.sequence.is_empty() {
    slot.node.seq.sequence = root_sequence.clone();
  }

  Ok(())
}

fn compute_msg_to_child(
  child: &SparseNodePartition,
  parent: &SparseNodePartition,
  edge: &SparseEdgePartition,
) -> Result<SparseSeqDistribution, Report> {
  let mut seq_dis = SparseSeqDistribution {
    variable: btreemap! {},
    variable_indel: BTreeSet::new(),
    fixed: btreemap! {},
    fixed_counts: parent.seq.composition.clone(),
    log_lh: forward_log_lh_remove_child(parent.profile.log_lh, edge.msg_from_child.log_lh),
  };
  let child_dis = &edge.msg_from_child;
  let mut parent_states = btreemap! {};
  let mut child_states = btreemap! {};
  for mutation in edge.fitch_subs() {
    child_states.insert(mutation.pos(), mutation.qry());
    parent_states.insert(mutation.pos(), mutation.reff());
  }
  for (pos, profile) in &parent.profile.variable {
    if !range_contains(&child.seq.non_char, *pos) {
      child_states.entry(*pos).or_insert(profile.state);
      parent_states.entry(*pos).or_insert(profile.state);
    }
  }
  for (pos, profile) in &child_dis.variable {
    if !range_contains(&child.seq.non_char, *pos) {
      child_states.entry(*pos).or_insert(profile.state);
      parent_states.entry(*pos).or_insert(profile.state);
    }
  }

  let mut delta_ll = LogLh::ZERO;
  for (pos, parent_state) in parent_states {
    let divisor = child_dis
      .variable
      .get(&pos)
      .map_or(&child_dis.fixed[&child_states[&pos]], |distribution| &distribution.dis);
    let numerator = parent
      .profile
      .variable
      .get(&pos)
      .map_or(&parent.profile.fixed[&parent_state], |distribution| &distribution.dis);
    let safe_divisor = divisor.mapv(|value| value.max(f64::MIN_POSITIVE));
    let mut dis = numerator / &safe_divisor;
    let normalization = normalize_1d_inplace(&mut dis, 1.0);
    delta_ll = forward_log_lh_add_normalization(delta_ll, normalization);
    seq_dis.variable.insert(
      pos,
      VarPos {
        dis,
        state: parent_state,
      },
    );
    seq_dis.fixed_counts.adjust_count(parent_state, -1);
  }
  for (state, profile) in &parent.profile.fixed {
    let safe_fixed = child_dis.fixed[state].mapv(|value| value.max(f64::MIN_POSITIVE));
    let mut dis = profile / &safe_fixed;
    let weight = seq_dis.fixed_counts.get(*state).unwrap() as f64;
    let normalization = normalize_1d_inplace(&mut dis, weight);
    delta_ll = forward_log_lh_add_normalization(delta_ll, normalization);
    seq_dis.fixed.insert(*state, dis);
  }
  seq_dis.log_lh += delta_ll;
  Ok(seq_dis)
}

fn compute_ml_subs_for_nodes(
  alphabet: &crate::alphabet::alphabet::Alphabet,
  parent: &SparseNodePartition,
  child: &SparseNodePartition,
  edge: &SparseEdgePartition,
) -> Result<Vec<Sub>, Report> {
  let positions = edge
    .fitch_subs()
    .iter()
    .map(Sub::pos)
    .chain(parent.profile.variable.keys().copied())
    .chain(child.profile.variable.keys().copied())
    .sorted()
    .dedup();
  positions
    .filter_map(|pos| {
      // Parity with the dense `edge_subs()` (`marginal_dense.rs`), which skips any position that
      // is `non_char` at either endpoint. A deleted position can still hold a `profile.variable`
      // entry whose argmax is an ordinary residue, and reporting it would put a substitution and
      // a deletion on the same edge at the same site.
      if range_contains(&parent.seq.non_char, pos) || range_contains(&child.seq.non_char, pos) {
        return None;
      }
      let parent_state = resolve_map_state(parent, pos, alphabet);
      let child_state = resolve_map_state(child, pos, alphabet);
      (parent_state != child_state && alphabet.is_canonical(parent_state) && alphabet.is_canonical(child_state))
        .then(|| Sub::new(parent_state, pos, child_state))
    })
    .collect()
}

fn resolve_map_state(
  node: &SparseNodePartition,
  pos: usize,
  alphabet: &crate::alphabet::alphabet::Alphabet,
) -> AsciiChar {
  if let Some(var) = node.profile.variable.get(&pos) {
    alphabet.char(argmax_first(&var.dis.view()).unwrap_or(0))
  } else {
    node.seq.sequence.get(pos).copied().unwrap_or_else(|| alphabet.char(0))
  }
}
