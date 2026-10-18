use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::partition::marginal::shared::update::MarginalBackward;
use crate::partition::marginal::sparse::candidate_sites::{CandidateSiteEvents, SiteSlot};
use crate::partition::marginal::sparse::message::{combine_messages, propagate_raw, propagate_raw_per_site};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::sparse::{
  FixedProfiles, SparseEdgeBackward, SparseEdgeObs, SparseNodeObs, SparseNodeState, SparseSeqDistribution,
};
use crate::partition::storage::var_pos_map::VarPosMap;
use eyre::Report;
use itertools::izip;
use ndarray::aview_mut1;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{GraphPass, GraphPassBackwardContext, GraphPassNodeOutput};
use treetime_primitives::LogLh;
use treetime_utils::interval::range::range_contains;

pub(crate) fn process_backward_indexed(
  partition: &PartitionMarginalSparse,
  gtr: &GTR,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
) -> Result<MarginalBackward<SparseNodeState, SparseEdgeBackward>, Report> {
  let pass = GraphPass::new(graph)?;
  let outputs = pass.map_backward(
    node_states,
    &partition.obs_edges,
    |key| treetime_utils::make_internal_error!("Partition node {key} is missing before the sparse marginal pass"),
    |context| process_node_backward_indexed(partition, gtr, branch_lengths, &context),
  )?;
  Ok(MarginalBackward {
    node_states: outputs.nodes,
    backward: outputs.edges,
  })
}

fn process_node_backward_indexed(
  partition: &PartitionMarginalSparse,
  gtr: &GTR,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  context: &GraphPassBackwardContext<'_, &SparseNodeState, SparseEdgeObs, SparseNodeState, SparseEdgeBackward>,
) -> Result<GraphPassNodeOutput<SparseNodeState, SparseEdgeBackward>, Report> {
  let alphabet = &partition.alphabet;
  let total_length = partition.total_length;
  let obs = &partition.obs_nodes[&context.key];

  let msg_to_parent = if context.is_leaf {
    leaf_message(obs, alphabet)?
  } else {
    combine_child_messages(partition, gtr, context, obs)?
  };

  let sequence = Arc::clone(&context.input.sequence);
  let (profile, parent_message) = if let Some((edge_key, edge_obs)) = context.parent_edge {
    let branch_length = fix_branch_length(total_length, branch_lengths[&edge_key]);
    let msg_from_child = if let Some(site_rates) = gtr.site_rates() {
      propagate_raw_per_site(
        gtr,
        site_rates,
        branch_length,
        true,
        &msg_to_parent,
        edge_obs.transmission.as_deref(),
      )
    } else {
      propagate_raw(
        &gtr.expQt(branch_length).t().to_owned(),
        &msg_to_parent,
        edge_obs.transmission.as_deref(),
      )
    };
    (
      SparseSeqDistribution::default(),
      Some(SparseEdgeBackward {
        msg_to_parent,
        msg_from_child,
      }),
    )
  } else {
    (msg_to_parent, None)
  };

  Ok(GraphPassNodeOutput {
    node: SparseNodeState { sequence, profile },
    parent_message,
  })
}

fn leaf_message(obs: &SparseNodeObs, alphabet: &Alphabet) -> Result<SparseSeqDistribution, Report> {
  let mut fixed = FixedProfiles::default();
  for state in alphabet.determined() {
    let count = obs.composition.get(state).unwrap_or(0);
    fixed.insert(state, count, alphabet.get_profile(state)?.iter().copied().collect());
  }
  let mut variable = VarPosMap::with_capacity(obs.fitch.variable.len());
  for (pos, profile) in &obs.fitch.variable {
    let dis = alphabet.construct_profile(profile.chars())?;
    let state = obs
      .fitch
      .chosen_state
      .get(pos)
      .copied()
      .filter(|state| alphabet.is_canonical(*state))
      .unwrap_or_else(|| profile.get_one());
    variable.insert_with(*pos, state, dis.len(), |row| aview_mut1(row).assign(&dis));
  }
  Ok(SparseSeqDistribution {
    variable,
    fixed,
    log_lh: LogLh::ZERO,
  })
}

fn combine_child_messages(
  partition: &PartitionMarginalSparse,
  gtr: &GTR,
  context: &GraphPassBackwardContext<'_, &SparseNodeState, SparseEdgeObs, SparseNodeState, SparseEdgeBackward>,
  obs: &SparseNodeObs,
) -> Result<SparseSeqDistribution, Report> {
  let alphabet = &partition.alphabet;
  let n_children = context.children.len();
  let mut child_messages = Vec::with_capacity(n_children);
  let mut n_events = 0;
  for child in context.children {
    let child_backward = child.edge_message()?;
    n_events +=
      2 * partition.obs_edges[&child.edge_key].fitch_subs().len() + child_backward.msg_from_child.variable.len();
    child_messages.push(&child_backward.msg_from_child);
  }

  let mut events = CandidateSiteEvents::with_capacity(n_children, n_events);
  for (ci, (child, child_message)) in izip!(context.children, &child_messages).enumerate() {
    for mutation in partition.obs_edges[&child.edge_key].fitch_subs() {
      events.insert(mutation.pos(), SiteSlot::State, mutation.reff());
      events.insert(mutation.pos(), SiteSlot::Reference(ci), mutation.qry());
    }
    for (pos, profile) in &child_message.variable {
      events.or_insert(*pos, SiteSlot::State, profile.state);
    }
  }
  let mut sites = events.finish();
  sites.fill_missing_references(|ci, pos, parent_state| {
    let child_obs = &partition.obs_nodes[&context.children[ci].node_key];
    if range_contains(&child_obs.non_char, pos) {
      if range_contains(&child_obs.gaps, pos) {
        alphabet.gap()
      } else {
        alphabet.unknown()
      }
    } else {
      parent_state
    }
  });

  combine_messages(
    &obs.composition,
    &child_messages,
    &sites,
    alphabet,
    context.is_root.then_some(gtr.pi()),
  )
}
