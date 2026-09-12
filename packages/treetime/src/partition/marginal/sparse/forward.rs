use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::partition::marginal::shared::normalize::{forward_log_lh_add_normalization, forward_log_lh_remove_child};
use crate::partition::marginal::sparse::message::{
  combine_messages, normalize_1d_inplace, propagate_raw, propagate_raw_per_site,
};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::marginal::sparse::reconstruct::{map_state, parsimony_seq};
use crate::partition::storage::sparse::{SparseEdgePartition, SparseNodePartition, SparseSeqDistribution, VarPos};
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::pass::{GraphPass, GraphPassForwardContext, GraphPassNodeOutput};
use treetime_primitives::{LogLh, Seq};
use treetime_utils::interval::range::range_contains;

pub fn process_forward_indexed(
  partition: &mut PartitionMarginalSparse,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<(), Report> {
  let alphabet = partition.alphabet.clone();
  let gtr = partition.gtr.clone();
  let length = partition.length;
  let root_sequence = partition.root_sequence.clone();
  let pass = GraphPass::new(graph)?;
  let outputs = pass.map_forward(
    &partition.nodes,
    &partition.edges,
    |key| treetime_utils::make_internal_error!("Partition node {key} is missing before the sparse marginal pass"),
    |context| process_node_forward_indexed(&alphabet, &gtr, length, &root_sequence, branch_lengths, &context),
  )?;
  partition.nodes = outputs.nodes;
  partition.edges = outputs.edges;
  Ok(())
}

fn process_node_forward_indexed(
  alphabet: &Alphabet,
  gtr: &GTR,
  length: usize,
  root_sequence: &Seq,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  context: &GraphPassForwardContext<'_, SparseNodePartition, SparseEdgePartition, SparseNodePartition>,
) -> Result<GraphPassNodeOutput<SparseNodePartition, SparseEdgePartition>, Report> {
  let mut node = context.input.clone();

  let parent_message = if let Some((edge_key, edge_data)) = context.parent_edge {
    let mut edge_data = edge_data.clone();
    let parent = context.parent.expect("Non-root node must have a parent");

    // Clone this node's parent-edge input and overwrite only the message and ML-subs fields. The edge
    // already carries `msg_from_child` from the backward pass and `fitch_subs`/`transmission` from the
    // Fitch pre-pass, and a fresh edge would corrupt the result.
    edge_data.msg_to_child = compute_msg_to_child(&node, parent, &edge_data)?;

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
      if !range_contains(&node.seq.non_char, *pos) {
        variable_pos.entry(*pos).or_insert(profile.state);
        parent_state.entry(*pos).or_insert(profile.state);
      }
    }
    for (pos, profile) in &edge_data.msg_to_parent.variable {
      variable_pos.entry(*pos).or_insert(profile.state);
      child_state.entry(*pos).or_insert(profile.state);
    }

    let branch_length = fix_branch_length(length, branch_lengths[&edge_key]);
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
    if context.is_leaf {
      edge_data.msg_from_parent = msg_from_parent.clone();
    }
    let profile = combine_messages(
      &node.seq.composition,
      &[msg_from_parent, edge_data.msg_to_parent.clone()],
      &variable_pos,
      &[parent_state, child_state],
      alphabet,
      None,
    )?;
    node.profile = profile;

    // Extend the parsimony chain. Leaves already hold their observed sequence, which is their
    // parsimony sequence; rebuilding it from the parent would discard the observed states a leaf
    // shares with its parent under Fitch compression.
    if !context.is_leaf && !parent.seq.sequence.is_empty() {
      node.seq.sequence = parsimony_seq(&parent.seq.sequence, &edge_data, &node, alphabet);
    }
    edge_data.set_ml_subs(compute_ml_subs_for_nodes(alphabet, parent, &node, &edge_data)?);

    Some(edge_data)
  } else {
    if node.seq.sequence.is_empty() {
      node.seq.sequence = root_sequence.clone();
    }
    None
  };

  Ok(GraphPassNodeOutput { node, parent_message })
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
  alphabet: &Alphabet,
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
      let parent_state = map_state(parent, pos, alphabet);
      let child_state = map_state(child, pos, alphabet);
      (parent_state != child_state && alphabet.is_canonical(parent_state) && alphabet.is_canonical(child_state))
        .then(|| Sub::new(parent_state, pos, child_state))
    })
    .collect()
}
