use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::partition::marginal::shared::normalize::{forward_log_lh_add_normalization, forward_log_lh_remove_child};
use crate::partition::marginal::shared::update::MarginalForward;
use crate::partition::marginal::sparse::message::{
  combine_messages, normalize_1d_inplace, propagate_raw, propagate_raw_per_site,
};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::marginal::sparse::reconstruct::{map_state, parsimony_seq};
use crate::partition::storage::sparse::{
  SparseEdgeBackward, SparseEdgeForward, SparseNodeObs, SparseNodeState, SparseSeqDistribution, VarPos,
};
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{GraphPass, GraphPassForwardContext, GraphPassNodeOutput};
use treetime_primitives::LogLh;
use treetime_utils::interval::range::range_contains;

/// Run the sparse marginal forward pass over borrowed inputs, node states, and backward messages,
/// returning updated node states, per-edge forward messages, and per-edge estimates (ML subs) as
/// distinct owned values.
pub fn process_forward_indexed(
  partition: &PartitionMarginalSparse,
  gtr: &GTR,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
  backward: &BTreeMap<GraphEdgeKey, SparseEdgeBackward>,
) -> Result<MarginalForward<SparseNodeState, SparseEdgeForward, Vec<Sub>>, Report> {
  let pass = GraphPass::new(graph)?;
  let outputs = pass.map_forward(
    node_states,
    backward,
    |key| treetime_utils::make_internal_error!("Partition node {key} is missing before the sparse marginal pass"),
    |context| process_node_forward_indexed(partition, gtr, branch_lengths, &context),
  )?;

  let mut forward = BTreeMap::new();
  let mut estimates = BTreeMap::new();
  for (edge_key, out) in outputs.edges {
    forward.insert(
      edge_key,
      SparseEdgeForward {
        msg_to_child: out.msg_to_child,
        msg_from_parent: out.msg_from_parent,
      },
    );
    estimates.insert(edge_key, out.subs_ml);
  }
  Ok(MarginalForward {
    node_states: outputs.nodes,
    forward,
    estimates,
  })
}

/// Combined per-edge output of the forward node visit, split by the driver into the distinct
/// forward-message and estimate owners.
struct SparseEdgeForwardOut {
  msg_to_child: SparseSeqDistribution,
  msg_from_parent: SparseSeqDistribution,
  subs_ml: Vec<Sub>,
}

fn process_node_forward_indexed(
  partition: &PartitionMarginalSparse,
  gtr: &GTR,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  context: &GraphPassForwardContext<'_, SparseNodeState, SparseEdgeBackward, SparseNodeState>,
) -> Result<GraphPassNodeOutput<SparseNodeState, SparseEdgeForwardOut>, Report> {
  let alphabet = &partition.alphabet;
  let length = partition.length;
  let obs = &partition.obs_nodes[&context.key];
  let mut node = context.input.clone();

  let parent_message = if let Some((edge_key, backward)) = context.parent_edge {
    let edge_obs = &partition.obs_edges[&edge_key];
    let parent = context.parent.expect("Non-root node must have a parent");
    let parent_key = context.parent_key.expect("Non-root node must have a parent key");
    let parent_obs = &partition.obs_nodes[&parent_key];

    let msg_to_child = compute_msg_to_child(obs, parent, parent_obs, edge_obs, backward)?;

    let mut variable_pos = btreemap! {};
    let mut parent_state = btreemap! {};
    let mut child_state = btreemap! {};
    for mutation in edge_obs.fitch_subs() {
      let current_state = mutation.qry();
      variable_pos.insert(mutation.pos(), current_state);
      parent_state.entry(mutation.pos()).or_insert_with(|| mutation.reff());
      child_state.entry(mutation.pos()).or_insert(current_state);
    }
    for (pos, profile) in &msg_to_child.variable {
      if !range_contains(&obs.non_char, *pos) {
        variable_pos.entry(*pos).or_insert(profile.state);
        parent_state.entry(*pos).or_insert(profile.state);
      }
    }
    for (pos, profile) in &backward.msg_to_parent.variable {
      variable_pos.entry(*pos).or_insert(profile.state);
      child_state.entry(*pos).or_insert(profile.state);
    }

    let branch_length = fix_branch_length(length, branch_lengths[&edge_key]);
    let msg_from_parent = if gtr.has_site_rates() {
      propagate_raw_per_site(
        gtr,
        branch_length,
        false,
        &msg_to_child,
        edge_obs.transmission.as_deref(),
      )
    } else {
      propagate_raw(
        &gtr.expQt(branch_length),
        &msg_to_child,
        edge_obs.transmission.as_deref(),
      )
    };
    // Persist the down-message only for tips: reconstruct_node_sequence imputes missing tip states
    // from it and has no branch length to recompute the propagation. Internal nodes never need it.
    let msg_from_parent_out = if context.is_leaf {
      msg_from_parent.clone()
    } else {
      SparseSeqDistribution::default()
    };
    let profile = combine_messages(
      &obs.composition,
      &[msg_from_parent, backward.msg_to_parent.clone()],
      &variable_pos,
      &[parent_state, child_state],
      alphabet,
      None,
    )?;
    node.profile = profile;

    // Extend the parsimony chain. Leaves already hold their observed sequence, which is their
    // parsimony sequence; rebuilding it from the parent would discard the observed states a leaf
    // shares with its parent under Fitch compression.
    if !context.is_leaf && !parent.sequence.is_empty() {
      node.sequence = parsimony_seq(&parent.sequence, edge_obs, obs, alphabet);
    }
    let subs_ml = compute_ml_subs_for_nodes(alphabet, parent, parent_obs, &node, obs, edge_obs)?;

    Some(SparseEdgeForwardOut {
      msg_to_child,
      msg_from_parent: msg_from_parent_out,
      subs_ml,
    })
  } else {
    if node.sequence.is_empty() {
      node.sequence = partition.root_sequence.clone();
    }
    None
  };

  Ok(GraphPassNodeOutput { node, parent_message })
}

fn compute_msg_to_child(
  child_obs: &SparseNodeObs,
  parent: &SparseNodeState,
  parent_obs: &SparseNodeObs,
  edge_obs: &crate::partition::storage::sparse::SparseEdgeObs,
  backward: &SparseEdgeBackward,
) -> Result<SparseSeqDistribution, Report> {
  let mut seq_dis = SparseSeqDistribution {
    variable: btreemap! {},
    variable_indel: BTreeSet::new(),
    fixed: btreemap! {},
    fixed_counts: parent_obs.composition.clone(),
    log_lh: forward_log_lh_remove_child(parent.profile.log_lh, backward.msg_from_child.log_lh),
  };
  let child_dis = &backward.msg_from_child;
  let mut parent_states = btreemap! {};
  let mut child_states = btreemap! {};
  for mutation in edge_obs.fitch_subs() {
    child_states.insert(mutation.pos(), mutation.qry());
    parent_states.insert(mutation.pos(), mutation.reff());
  }
  for (pos, profile) in &parent.profile.variable {
    if !range_contains(&child_obs.non_char, *pos) {
      child_states.entry(*pos).or_insert(profile.state);
      parent_states.entry(*pos).or_insert(profile.state);
    }
  }
  for (pos, profile) in &child_dis.variable {
    if !range_contains(&child_obs.non_char, *pos) {
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
  parent: &SparseNodeState,
  parent_obs: &SparseNodeObs,
  child: &SparseNodeState,
  child_obs: &SparseNodeObs,
  edge_obs: &crate::partition::storage::sparse::SparseEdgeObs,
) -> Result<Vec<Sub>, Report> {
  let positions = edge_obs
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
      if range_contains(&parent_obs.non_char, pos) || range_contains(&child_obs.non_char, pos) {
        return None;
      }
      let parent_state = map_state(parent, pos, alphabet);
      let child_state = map_state(child, pos, alphabet);
      (parent_state != child_state && alphabet.is_canonical(parent_state) && alphabet.is_canonical(child_state))
        .then(|| Sub::new(parent_state, pos, child_state))
    })
    .collect()
}
