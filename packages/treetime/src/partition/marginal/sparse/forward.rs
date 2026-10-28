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
use crate::partition::marginal::sparse::reconstruct::{changed_group_positions, map_state};
use crate::partition::marginal::sparse::site_grouping::{GroupedSites, ProfileSource, group_layer, group_sites};
use crate::partition::storage::site_groups::SiteGroups;
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
use treetime_utils::interval::range::range_contains_sorted;
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
    } = compute_msg_to_child(
      alphabet,
      &node.sequence,
      parent,
      parent_obs,
      edge_obs,
      backward,
      partition.impute_missing_data,
    )?;

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
      if !range_contains_sorted(&backward.child_no_data, *pos) {
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
    let msg_from_parent = propagate_forward(gtr, branch_length, &msg_to_child)?;
    let groups = group_sites(
      &[
        group_layer(&msg_from_parent.groups),
        group_layer(&backward.msg_to_parent.groups),
      ],
      &[],
      &sites,
      &node.sequence,
      alphabet,
    );
    let profile = combine_messages(
      &obs.composition,
      &[&msg_from_parent, &backward.msg_to_parent],
      &sites,
      groups,
      alphabet,
      None,
    )?;
    node.profile = profile;
    let mut msg_from_parent = if is_leaf {
      msg_from_parent
    } else {
      SparseSeqDistribution::default()
    };
    if !child_missing.is_empty() {
      let child_missing = propagate_forward(
        gtr,
        branch_length,
        &SparseSeqDistribution {
          variable: child_missing,
          ..SparseSeqDistribution::default()
        },
      )?;
      let outside = if is_leaf {
        &mut msg_from_parent.variable
      } else {
        &mut node.profile.variable
      };
      for (&pos, entry) in &child_missing.variable {
        outside.insert(pos, entry.state, entry.dis);
      }
    }

    let subs_ml = compute_ml_subs_for_nodes(alphabet, parent, parent_obs, node, obs, edge_obs, backward)?;

    Ok(SparseEdgeForwardOut {
      msg_to_child,
      msg_from_parent,
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

fn propagate_forward(
  gtr: &GTR,
  branch_length: f64,
  message: &SparseSeqDistribution,
) -> Result<SparseSeqDistribution, Report> {
  if let Some(site_rates) = gtr.site_rates() {
    propagate_raw_per_site(gtr, site_rates, branch_length, false, message)
  } else {
    Ok(propagate_raw(&gtr.expQt(branch_length), message))
  }
}

fn compute_msg_to_child(
  alphabet: &Alphabet,
  child_sequence: &Seq,
  parent: &SparseNodeState,
  parent_obs: &SparseNodeObs,
  edge_obs: &SparseEdgeObs,
  backward: &SparseEdgeBackward,
  impute_missing_data: bool,
) -> Result<OutsideMessage, Report> {
  let child_dis = &backward.msg_from_child;
  let n_events = 2 * (edge_obs.fitch_subs().len() + parent.profile.variable.len() + child_dis.variable.len());
  let mut events = CandidateSiteEvents::with_capacity(1, n_events);
  for mutation in edge_obs.fitch_subs() {
    events.insert(mutation.pos(), SiteSlot::Reference(0), mutation.qry());
    events.insert(mutation.pos(), SiteSlot::State, mutation.reff());
  }
  for (pos, profile) in parent.profile.variable.iter().chain(&child_dis.variable) {
    if !range_contains_sorted(&backward.child_no_data, *pos) {
      events.or_insert(*pos, SiteSlot::Reference(0), profile.state);
      events.or_insert(*pos, SiteSlot::State, profile.state);
    }
  }
  let mut sites = events.finish();
  sites.relabel(|pos, state| canonical_or(alphabet, parent.sequence[pos], state));

  let mut division = PosteriorDivision::new(alphabet, parent, child_dis, sites.len());
  for (pos, parent_state, child_state) in sites.iter() {
    let child_state = child_state[0].ok_or_else(|| make_internal_report!("Position {pos} has no child state"))?;
    division.divide_site(
      pos,
      parent_state,
      canonical_or(alphabet, child_sequence[pos], child_state),
    )?;
  }
  let groups = group_sites(
    &[group_layer(&parent.profile.groups), group_layer(&child_dis.groups)],
    &backward.child_no_data,
    &sites,
    &parent.sequence,
    alphabet,
  );
  let group_profiles = groups
    .into_iter()
    .map(|group| division.divide_group(group))
    .filter_ok(|(_, profiles)| !profiles.is_empty())
    .collect::<Result<Vec<_>, Report>>()?;
  division.msg_to_child.groups = SiteGroups::from_groups(group_profiles);
  let child_missing = division.take_child_missing(
    sites.non_candidate_positions(&backward.child_no_data),
    impute_missing_data,
  )?;
  division.divide_fixed(parent_obs)?;
  Ok(OutsideMessage {
    msg_to_child: division.finish(),
    child_missing,
  })
}

struct PosteriorDivision<'a> {
  alphabet: &'a Alphabet,
  parent: &'a SparseNodeState,
  child_dis: &'a SparseSeqDistribution,
  msg_to_child: SparseSeqDistribution,
  delta_ll: LogLh,
  removed: AsciiTable<usize>,
}

impl<'a> PosteriorDivision<'a> {
  fn new(
    alphabet: &'a Alphabet,
    parent: &'a SparseNodeState,
    child_dis: &'a SparseSeqDistribution,
    n_sites: usize,
  ) -> Self {
    let msg_to_child = SparseSeqDistribution {
      variable: VarPosMap::with_capacity(n_sites),
      groups: SiteGroups::default(),
      fixed: FixedProfiles::default(),
      log_lh: forward_log_lh_remove_child(parent.profile.log_lh, child_dis.log_lh),
    };
    Self {
      alphabet,
      parent,
      child_dis,
      msg_to_child,
      delta_ll: LogLh::ZERO,
      removed: AsciiTable::default(),
    }
  }

  fn divide_site(&mut self, pos: usize, parent_state: AsciiChar, child_state: AsciiChar) -> Result<(), Report> {
    let divisor = self
      .child_dis
      .profile_at(pos, child_state)
      .ok_or_else(|| make_internal_report!("Child message has no profile of state {child_state} at position {pos}"))?;
    let numerator = self.parent.profile.profile_at(pos, parent_state).ok_or_else(|| {
      make_internal_report!("Parent posterior has no profile of state {parent_state} at position {pos}")
    })?;
    let mut dis = divide_out_site(numerator, divisor);
    let normalization = normalize_1d_inplace(aview_mut1(&mut dis), 1.0);
    self.delta_ll = forward_log_lh_add_normalization(self.delta_ll, normalization);
    if self.alphabet.is_canonical(self.parent.sequence[pos]) {
      self.msg_to_child.variable.insert(pos, parent_state, &dis);
      self.removed[parent_state] += 1;
    } else {
      self.msg_to_child.variable.insert(pos, child_state, &dis);
    }
    Ok(())
  }

  fn divide_group(
    &mut self,
    GroupedSites {
      sources,
      ranges,
      counts,
    }: GroupedSites,
  ) -> Result<(Vec<(usize, usize)>, FixedProfiles), Report> {
    let mut profiles = FixedProfiles::default();
    for (state, count) in counts {
      let numerator = group_entry(&self.parent.profile, sources[0], state);
      let divisor = group_entry(self.child_dis, sources[1], state);
      if numerator.is_none() && divisor.is_none() {
        continue;
      }
      let numerator = shared_or_fixed(numerator, &self.parent.profile, state)?;
      let divisor = shared_or_fixed(divisor, self.child_dis, state)?;
      let mut dis = divide_out_site(numerator, divisor);
      let normalization = normalize_1d_inplace(aview_mut1(&mut dis), count.to_f64_exact());
      self.delta_ll = forward_log_lh_add_normalization(self.delta_ll, normalization);
      self.removed[state] += count;
      profiles.insert(state, count, dis);
    }
    Ok((ranges, profiles))
  }

  fn take_child_missing(
    &mut self,
    positions: impl Iterator<Item = usize>,
    impute_missing_data: bool,
  ) -> Result<VarPosMap, Report> {
    let mut child_missing = VarPosMap::default();
    for pos in positions {
      let parent_state = self.parent.sequence[pos];
      if self.alphabet.is_canonical(parent_state) {
        self.removed[parent_state] += 1;
        if impute_missing_data {
          let posterior = self.parent.profile.profile_at(pos, parent_state).ok_or_else(|| {
            make_internal_report!("Parent posterior has no profile of state {parent_state} at position {pos}")
          })?;
          child_missing.insert(pos, parent_state, posterior);
        }
      } else if impute_missing_data && let Some(posterior) = self.parent.profile.variable.get(pos) {
        child_missing.insert(pos, posterior.state, posterior.dis);
      }
    }
    Ok(child_missing)
  }

  fn divide_fixed(&mut self, parent_obs: &SparseNodeObs) -> Result<(), Report> {
    for fixed in self.parent.profile.fixed.iter() {
      let state = fixed.state;
      let mut dis = divide_out_site(&fixed.dis, &self.child_dis.fixed[state].dis);
      let count = parent_obs
        .composition
        .get(state)
        .ok_or_else(|| make_internal_report!("Fixed state {state} has no count"))?
        .saturating_sub(self.removed[state]);
      let normalization = normalize_1d_inplace(aview_mut1(&mut dis), count.to_f64_exact());
      self.delta_ll = forward_log_lh_add_normalization(self.delta_ll, normalization);
      self.msg_to_child.fixed.insert(state, count, dis);
    }
    Ok(())
  }

  fn finish(mut self) -> SparseSeqDistribution {
    self.msg_to_child.log_lh += self.delta_ll;
    self.msg_to_child
  }
}

fn group_entry(message: &SparseSeqDistribution, source: ProfileSource, state: AsciiChar) -> Option<&[f64]> {
  match source {
    ProfileSource::Group(group) => message
      .groups
      .profiles(group)
      .get(state)
      .map(|entry| entry.dis.as_slice()),
    ProfileSource::Fixed | ProfileSource::Missing => None,
  }
}

fn shared_or_fixed<'a>(
  entry: Option<&'a [f64]>,
  message: &'a SparseSeqDistribution,
  state: AsciiChar,
) -> Result<&'a [f64], Report> {
  entry
    .or_else(|| message.fixed.get(state).map(|fixed| fixed.dis.as_slice()))
    .ok_or_else(|| make_internal_report!("Message has no fixed profile of state {state}"))
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
  backward: &SparseEdgeBackward,
) -> Result<Vec<Sub>, Report> {
  let positions = edge_obs
    .fitch_subs()
    .iter()
    .map(Sub::pos)
    .chain(parent.profile.variable.keys().copied())
    .chain(child.profile.variable.keys().copied())
    .chain(changed_group_positions(parent, alphabet))
    .chain(changed_group_positions(child, alphabet))
    .sorted_unstable()
    .dedup();
  positions
    .filter_map(|pos| {
      if range_contains_sorted(&parent_obs.non_char, pos)
        || range_contains_sorted(&child_obs.non_char, pos)
        || range_contains_sorted(&backward.child_no_data, pos)
      {
        return None;
      }
      let parent_state = map_state(parent, pos, alphabet);
      let child_state = map_state(child, pos, alphabet);
      (parent_state != child_state && alphabet.is_canonical(parent_state) && alphabet.is_canonical(child_state))
        .then(|| Sub::new(parent_state, pos, child_state))
    })
    .collect()
}
