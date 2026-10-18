use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::partition::marginal::shared::normalize::{forward_log_lh_add_normalization, forward_log_lh_remove_child};
use crate::partition::marginal::shared::update::MarginalForward;
use crate::partition::marginal::sparse::candidate_sites::{CandidateSiteEvents, SiteSlot};
use crate::partition::marginal::sparse::message::{
  combine_messages, divide_out_site, normalize_1d_inplace, propagate_raw, propagate_raw_per_site,
};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::marginal::sparse::reconstruct::{map_state, parsimony_seq};
use crate::partition::storage::sparse::{
  FixedProfiles, SparseEdgeBackward, SparseEdgeForward, SparseEdgeObs, SparseNodeObs, SparseNodeState,
  SparseSeqDistribution,
};
use crate::partition::storage::var_pos_map::VarPosMap;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use ndarray::aview_mut1;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{GraphPass, GraphPassForwardContext, GraphPassNodeOutput};
use treetime_primitives::LogLh;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::interval::range::range_contains;
use treetime_utils::make_internal_report;

pub(crate) fn process_forward_indexed(
  partition: &PartitionMarginalSparse,
  gtr: &GTR,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  node_states: BTreeMap<GraphNodeKey, SparseNodeState>,
  backward: &BTreeMap<GraphEdgeKey, SparseEdgeBackward>,
) -> Result<MarginalForward<SparseNodeState, SparseEdgeForward, Vec<Sub>>, Report> {
  let pass = GraphPass::new(graph)?;
  let outputs = pass.map_forward_owned(
    node_states,
    backward,
    |key| treetime_utils::make_internal_error!("Partition node {key} is missing before the sparse marginal pass"),
    |context| process_node_forward_indexed(partition, gtr, branch_lengths, context),
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

fn process_node_forward_indexed(
  partition: &PartitionMarginalSparse,
  gtr: &GTR,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  context: GraphPassForwardContext<'_, SparseNodeState, SparseEdgeBackward, SparseNodeState>,
) -> Result<GraphPassNodeOutput<SparseNodeState, SparseEdgeForwardOut>, Report> {
  let alphabet = &partition.alphabet;
  let total_length = partition.total_length;
  let obs = &partition.obs_nodes[&context.key];
  let mut node = context.input;

  let parent_message = if let (Some((edge_key, backward)), Some(parent), Some(parent_key)) =
    (context.parent_edge, context.parent, context.parent_key)
  {
    let edge_obs = &partition.obs_edges[&edge_key];
    let parent_obs = &partition.obs_nodes[&parent_key];

    let msg_to_child = compute_msg_to_child(obs, parent, parent_obs, edge_obs, backward)?;

    let n_events =
      3 * edge_obs.fitch_subs().len() + 2 * (msg_to_child.variable.len() + backward.msg_to_parent.variable.len());
    let mut events = CandidateSiteEvents::with_capacity(2, n_events);
    for mutation in edge_obs.fitch_subs() {
      let current_state = mutation.qry();
      events.insert(mutation.pos(), SiteSlot::State, current_state);
      events.or_insert(mutation.pos(), SiteSlot::Reference(0), mutation.reff());
      events.or_insert(mutation.pos(), SiteSlot::Reference(1), current_state);
    }
    for (pos, profile) in &msg_to_child.variable {
      if !range_contains(&obs.non_char, *pos) {
        events.or_insert(*pos, SiteSlot::State, profile.state);
        events.or_insert(*pos, SiteSlot::Reference(0), profile.state);
      }
    }
    for (pos, profile) in &backward.msg_to_parent.variable {
      events.or_insert(*pos, SiteSlot::State, profile.state);
      events.or_insert(*pos, SiteSlot::Reference(1), profile.state);
    }
    let sites = events.finish();

    let branch_length = fix_branch_length(total_length, branch_lengths[&edge_key]);
    let msg_from_parent = if let Some(site_rates) = gtr.site_rates() {
      propagate_raw_per_site(
        gtr,
        site_rates,
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
    let profile = combine_messages(
      &obs.composition,
      &[&msg_from_parent, &backward.msg_to_parent],
      &sites,
      alphabet,
      None,
    )?;
    node.profile = profile;
    let msg_from_parent_out = if context.is_leaf {
      msg_from_parent
    } else {
      SparseSeqDistribution::default()
    };

    if !context.is_leaf && !parent.sequence.is_empty() {
      let sequence = parsimony_seq(&parent.sequence, edge_obs, obs, alphabet);
      if *node.sequence != sequence {
        node.sequence = Arc::new(sequence);
      }
    }
    let subs_ml = compute_ml_subs_for_nodes(alphabet, parent, parent_obs, &node, obs, edge_obs)?;

    Some(SparseEdgeForwardOut {
      msg_to_child,
      msg_from_parent: msg_from_parent_out,
      subs_ml,
    })
  } else {
    if node.sequence.is_empty() {
      node.sequence = Arc::new(partition.root_sequence.clone());
    }
    None
  };

  Ok(GraphPassNodeOutput { node, parent_message })
}

struct SparseEdgeForwardOut {
  msg_to_child: SparseSeqDistribution,
  msg_from_parent: SparseSeqDistribution,
  subs_ml: Vec<Sub>,
}

fn compute_msg_to_child(
  child_obs: &SparseNodeObs,
  parent: &SparseNodeState,
  parent_obs: &SparseNodeObs,
  edge_obs: &SparseEdgeObs,
  backward: &SparseEdgeBackward,
) -> Result<SparseSeqDistribution, Report> {
  let child_dis = &backward.msg_from_child;
  let n_events = 2 * (edge_obs.fitch_subs().len() + parent.profile.variable.len() + child_dis.variable.len());
  let mut events = CandidateSiteEvents::with_capacity(1, n_events);
  for mutation in edge_obs.fitch_subs() {
    events.insert(mutation.pos(), SiteSlot::Reference(0), mutation.qry());
    events.insert(mutation.pos(), SiteSlot::State, mutation.reff());
  }
  for (pos, profile) in parent.profile.variable.iter().chain(&child_dis.variable) {
    if !range_contains(&child_obs.non_char, *pos) {
      events.or_insert(*pos, SiteSlot::Reference(0), profile.state);
      events.or_insert(*pos, SiteSlot::State, profile.state);
    }
  }
  let sites = events.finish();

  let mut seq_dis = SparseSeqDistribution {
    variable: VarPosMap::with_capacity(sites.len()),
    fixed: FixedProfiles::default(),
    log_lh: forward_log_lh_remove_child(parent.profile.log_lh, backward.msg_from_child.log_lh),
  };
  let mut delta_ll = LogLh::ZERO;
  let mut removed = [0_usize; 128];
  for (pos, parent_state, child_state) in sites.iter() {
    let child_state = child_state[0].ok_or_else(|| make_internal_report!("Position {pos} has no child state"))?;
    let divisor = child_dis.variable.get(pos).map_or_else(
      || child_dis.fixed[child_state].dis.as_slice(),
      |distribution| distribution.dis,
    );
    let numerator = parent.profile.variable.get(pos).map_or_else(
      || parent.profile.fixed[parent_state].dis.as_slice(),
      |distribution| distribution.dis,
    );
    let mut dis = divide_out_site(numerator, divisor);
    let normalization = normalize_1d_inplace(aview_mut1(&mut dis), 1.0);
    delta_ll = forward_log_lh_add_normalization(delta_ll, normalization);
    seq_dis.variable.insert(pos, parent_state, &dis);
    removed[usize::from(parent_state)] += 1;
  }
  for fixed in parent.profile.fixed.iter() {
    let state = fixed.state;
    let mut dis = divide_out_site(&fixed.dis, &child_dis.fixed[state].dis);
    let count = parent_obs
      .composition
      .get(state)
      .ok_or_else(|| make_internal_report!("Fixed state {state} has no count"))?
      .saturating_sub(removed[usize::from(state)]);
    let normalization = normalize_1d_inplace(aview_mut1(&mut dis), count.to_f64_exact());
    delta_ll = forward_log_lh_add_normalization(delta_ll, normalization);
    seq_dis.fixed.insert(state, count, dis);
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
  edge_obs: &SparseEdgeObs,
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
