use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::partition::marginal::shared::normalize::{forward_log_lh_add_normalization, forward_log_lh_remove_child};
use crate::partition::marginal::shared::update::MarginalForward;
use crate::partition::marginal::sparse::backward::own_state;
use crate::partition::marginal::sparse::candidate_sites::{CandidateSiteEvents, SiteSlot};
use crate::partition::marginal::sparse::message::{
  combine_messages, divide_out_site, normalize_1d_inplace, propagate_raw, propagate_raw_per_site,
};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::marginal::sparse::reconstruct::map_state;
use crate::partition::storage::sparse::{
  FixedProfiles, SparseEdgeBackward, SparseEdgeForward, SparseEdgeObs, SparseNodeObs, SparseNodeState,
  SparseSeqDistribution,
};
use crate::partition::storage::var_pos_map::VarPosMap;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use ndarray::aview_mut1;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::schedule_scan::{NodeEdgeMaps, PerEdge};
use treetime_primitives::{AsciiChar, AsciiTable, LogLh, Seq};
use treetime_utils::cast::ToF64Exact;
use treetime_utils::interval::range::range_contains;
use treetime_utils::make_internal_report;

pub(crate) fn process_forward_indexed(
  partition: &PartitionMarginalSparse,
  gtr: &GTR,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  node_states: NodeMap<SparseNodeState>,
  backward: &EdgeMap<SparseEdgeBackward>,
) -> Result<MarginalForward<SparseNodeState, SparseEdgeForward, Vec<Sub>>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  let pass = SparseForwardPass {
    schedule: &schedule,
    partition,
    gtr,
    branch_lengths,
  };
  let backward = backward.iter().collect();
  let outputs = schedule.scan_down(
    schedule.join_parent_edges(node_states, backward)?,
    |key, (node, mut backward), parents: PerEdge<SparseForwardOutput>| {
      pass.node(key, node, backward.pop(), parents.first().map(|parent| &parent.0))
    },
  )?;
  let NodeEdgeMaps { nodes, edges } = outputs.unwrap_shared().split_parent_edges()?;

  let mut forward = EdgeMap::new();
  let mut estimates = EdgeMap::new();
  for (edge_key, out) in edges {
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
    node_states: nodes,
    forward,
    estimates,
  })
}

struct SparseForwardPass<'a> {
  schedule: &'a GraphSchedule<'a>,
  partition: &'a PartitionMarginalSparse,
  gtr: &'a GTR,
  branch_lengths: &'a EdgeMap<f64>,
}

impl SparseForwardPass<'_> {
  fn node(
    &self,
    key: GraphNodeKey,
    mut node: SparseNodeState,
    backward: Option<&SparseEdgeBackward>,
    parent: Option<&SparseNodeState>,
  ) -> Result<(SparseForwardOutput, SparseForwardOutput), Report> {
    let partition = self.partition;
    let edge = match (self.schedule.parents(key).first(), backward, parent) {
      (Some(&(parent_key, edge_key)), Some(backward), Some(parent)) => {
        Some(self.edge(key, parent_key, edge_key, &mut node, backward, parent)?)
      },
      (None, None, None) => None,
      _ => {
        return Err(make_internal_report!(
          "Sparse marginal forward pass received an inconsistent parent edge and parent message at node {key}"
        ));
      },
    };
    let output = Arc::new((node, edge.into_iter().collect()));
    Ok((Arc::clone(&output), output))
  }

  fn edge(
    &self,
    key: GraphNodeKey,
    parent_key: GraphNodeKey,
    edge_key: GraphEdgeKey,
    node: &mut SparseNodeState,
    backward: &SparseEdgeBackward,
    parent: &SparseNodeState,
  ) -> Result<SparseEdgeForwardOut, Report> {
    let partition = self.partition;
    let gtr = self.gtr;
    let alphabet = &partition.alphabet;
    let obs = &partition.obs_nodes[&key];
    let is_leaf = self.schedule.children(key).is_empty();
    let edge_obs = &partition.obs_edges[&edge_key];
    let parent_obs = &partition.obs_nodes[&parent_key];

    let OutsideMessage {
      msg_to_child,
      child_missing,
    } = compute_msg_to_child(alphabet, &node.sequence, parent, parent_obs, edge_obs, backward)?;

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
      if !range_contains(&backward.child_no_data, *pos) {
        events.or_insert(*pos, SiteSlot::State, profile.state);
        events.or_insert(*pos, SiteSlot::Reference(0), profile.state);
      }
    }
    for (pos, profile) in &backward.msg_to_parent.variable {
      events.or_insert(*pos, SiteSlot::State, profile.state);
      events.or_insert(*pos, SiteSlot::Reference(1), profile.state);
    }
    let mut sites = events.finish();
    sites.relabel(|pos, state| own_state(alphabet, obs, &node.sequence, pos, state));

    let branch_length = fix_branch_length(partition.total_length, self.branch_lengths[&edge_key]);
    let msg_from_parent = propagate_forward(gtr, branch_length, &msg_to_child);
    let profile = combine_messages(
      &obs.composition,
      &[&msg_from_parent, &backward.msg_to_parent],
      &sites,
      alphabet,
      None,
    )?;
    node.profile = profile;
    let msg_from_parent_out = if is_leaf {
      let child_missing = propagate_forward(
        gtr,
        branch_length,
        &SparseSeqDistribution {
          variable: child_missing,
          ..SparseSeqDistribution::default()
        },
      );
      let mut msg_from_parent = msg_from_parent;
      for (&pos, outside) in &child_missing.variable {
        msg_from_parent.variable.insert(pos, outside.state, outside.dis);
      }
      msg_from_parent
    } else {
      SparseSeqDistribution::default()
    };

    let subs_ml = compute_ml_subs_for_nodes(alphabet, parent, parent_obs, node, obs, edge_obs)?;

    Ok(SparseEdgeForwardOut {
      msg_to_child,
      msg_from_parent: msg_from_parent_out,
      subs_ml,
    })
  }
}

type SparseForwardOutput = Arc<(SparseNodeState, PerEdge<SparseEdgeForwardOut>)>;

struct SparseEdgeForwardOut {
  msg_to_child: SparseSeqDistribution,
  msg_from_parent: SparseSeqDistribution,
  subs_ml: Vec<Sub>,
}

fn propagate_forward(gtr: &GTR, branch_length: f64, message: &SparseSeqDistribution) -> SparseSeqDistribution {
  if let Some(site_rates) = gtr.site_rates() {
    propagate_raw_per_site(gtr, site_rates, branch_length, false, message)
  } else {
    propagate_raw(&gtr.expQt(branch_length), message)
  }
}

fn compute_msg_to_child(
  alphabet: &Alphabet,
  child_sequence: &Seq,
  parent: &SparseNodeState,
  parent_obs: &SparseNodeObs,
  edge_obs: &SparseEdgeObs,
  backward: &SparseEdgeBackward,
) -> Result<OutsideMessage, Report> {
  let child_dis = &backward.msg_from_child;
  let n_events = 2 * (edge_obs.fitch_subs().len() + parent.profile.variable.len() + child_dis.variable.len());
  let mut events = CandidateSiteEvents::with_capacity(1, n_events);
  for mutation in edge_obs.fitch_subs() {
    events.insert(mutation.pos(), SiteSlot::Reference(0), mutation.qry());
    events.insert(mutation.pos(), SiteSlot::State, mutation.reff());
  }
  for (pos, profile) in parent.profile.variable.iter().chain(&child_dis.variable) {
    if !range_contains(&backward.child_no_data, *pos) {
      events.or_insert(*pos, SiteSlot::Reference(0), profile.state);
      events.or_insert(*pos, SiteSlot::State, profile.state);
    }
  }
  let mut sites = events.finish();
  sites.relabel(|pos, state| canonical_or(alphabet, parent.sequence[pos], state));

  let mut seq_dis = SparseSeqDistribution {
    variable: VarPosMap::with_capacity(sites.len()),
    fixed: FixedProfiles::default(),
    log_lh: forward_log_lh_remove_child(parent.profile.log_lh, backward.msg_from_child.log_lh),
  };
  let mut delta_ll = LogLh::ZERO;
  let mut removed = AsciiTable::<usize>::default();
  for (pos, parent_state, child_state) in sites.iter() {
    let child_state = child_state[0].ok_or_else(|| make_internal_report!("Position {pos} has no child state"))?;
    let child_state = canonical_or(alphabet, child_sequence[pos], child_state);
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
    if alphabet.is_canonical(parent.sequence[pos]) {
      seq_dis.variable.insert(pos, parent_state, &dis);
      removed[parent_state] += 1;
    } else {
      seq_dis.variable.insert(pos, child_state, &dis);
    }
  }
  let mut child_missing = VarPosMap::default();
  for pos in backward.child_no_data.iter().flat_map(|&(start, end)| start..end) {
    let parent_state = parent.sequence[pos];
    if sites.contains(pos) {
      continue;
    }
    if alphabet.is_canonical(parent_state) {
      let posterior = parent.profile.variable.get(pos).map_or_else(
        || parent.profile.fixed[parent_state].dis.as_slice(),
        |distribution| distribution.dis,
      );
      child_missing.insert(pos, parent_state, posterior);
      removed[parent_state] += 1;
    } else if let Some(posterior) = parent.profile.variable.get(pos) {
      child_missing.insert(pos, posterior.state, posterior.dis);
    }
  }
  for fixed in parent.profile.fixed.iter() {
    let state = fixed.state;
    let mut dis = divide_out_site(&fixed.dis, &child_dis.fixed[state].dis);
    let count = parent_obs
      .composition
      .get(state)
      .ok_or_else(|| make_internal_report!("Fixed state {state} has no count"))?
      .saturating_sub(removed[state]);
    let normalization = normalize_1d_inplace(aview_mut1(&mut dis), count.to_f64_exact());
    delta_ll = forward_log_lh_add_normalization(delta_ll, normalization);
    seq_dis.fixed.insert(state, count, dis);
  }
  seq_dis.log_lh += delta_ll;
  Ok(OutsideMessage {
    msg_to_child: seq_dis,
    child_missing,
  })
}

fn canonical_or(alphabet: &Alphabet, state: AsciiChar, fallback: AsciiChar) -> AsciiChar {
  if alphabet.is_canonical(state) { state } else { fallback }
}

struct OutsideMessage {
  msg_to_child: SparseSeqDistribution,
  child_missing: VarPosMap,
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
