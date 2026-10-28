use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::partition::marginal::shared::update::MarginalBackward;
use crate::partition::marginal::sparse::candidate_sites::{CandidateSiteEvents, SiteSlot};
use crate::partition::marginal::sparse::message::{combine_messages, propagate_raw, propagate_raw_per_site};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::marginal::sparse::sequences::fitch_node_sequences;
use crate::partition::storage::sparse::{
  FixedProfiles, SparseEdgeBackward, SparseEdgeObs, SparseNodeObs, SparseNodeState, SparseSeqDistribution,
};
use crate::partition::storage::var_pos_map::VarPosMap;
use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::aview_mut1;
use smallvec::smallvec;
use std::sync::Arc;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::schedule_scan::{NodeEdgeMaps, PerEdge};
use treetime_primitives::{LogLh, Seq};
use treetime_utils::interval::range::range_contains;
use treetime_utils::interval::range_difference::range_difference;
use treetime_utils::interval::range_intersection::range_intersection_iter;
use treetime_utils::interval::range_union::range_union;
use treetime_utils::make_internal_report;

pub(crate) fn process_backward_indexed(
  partition: &PartitionMarginalSparse,
  gtr: &GTR,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  node_states: &NodeMap<SparseNodeState>,
) -> Result<MarginalBackward<SparseNodeState, SparseEdgeBackward>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  let sequences = fitch_node_sequences(partition, &schedule, node_states)?;
  let pass = SparseBackwardPass {
    schedule: &schedule,
    partition,
    gtr,
    branch_lengths,
  };
  let edges = partition.obs_edges.iter().collect();
  let outputs = schedule.scan_up(
    schedule.join_parent_edges(sequences, edges)?,
    |key, (sequence, mut edge_obs), children| pass.node(key, sequence, edge_obs.pop(), &children),
  )?;
  let NodeEdgeMaps { nodes, edges } = outputs.unwrap_shared().split_parent_edges()?;
  Ok(MarginalBackward {
    node_states: nodes,
    backward: edges,
  })
}

struct SparseBackwardPass<'a> {
  schedule: &'a GraphSchedule<'a>,
  partition: &'a PartitionMarginalSparse,
  gtr: &'a GTR,
  branch_lengths: &'a EdgeMap<f64>,
}

impl SparseBackwardPass<'_> {
  fn node(
    &self,
    key: GraphNodeKey,
    sequence: Arc<Seq>,
    edge_obs: Option<&SparseEdgeObs>,
    children: &[SparseBackwardOutput],
  ) -> Result<(SparseBackwardOutput, SparseBackwardOutput), Report> {
    let partition = self.partition;
    let gtr = self.gtr;
    let obs = &partition.obs_nodes[&key];

    let (msg_to_parent, no_data) = if self.schedule.children(key).is_empty() {
      (
        leaf_message(obs, &partition.alphabet)?,
        range_union(&[obs.gaps.clone(), obs.unknown.clone()]),
      )
    } else {
      let child_edges = self.child_edges(key, children)?;
      let no_data = range_intersection_iter(child_edges.iter().map(|edge| &edge.child_no_data)).collect_vec();
      (self.combine_child_messages(key, &child_edges, obs, &sequence)?, no_data)
    };

    let (profile, parent_message) = match (self.schedule.parents(key).first(), edge_obs) {
      (Some(&(_, edge_key)), Some(_)) => {
        let branch_length = fix_branch_length(partition.total_length, self.branch_lengths[&edge_key]);
        let msg_from_child = if let Some(site_rates) = gtr.site_rates() {
          propagate_raw_per_site(gtr, site_rates, branch_length, true, &msg_to_parent)
        } else {
          propagate_raw(&gtr.expQt(branch_length).t().to_owned(), &msg_to_parent)
        };
        (
          SparseSeqDistribution::default(),
          smallvec![SparseEdgeBackward {
            msg_to_parent,
            msg_from_child,
            child_no_data: no_data,
          }],
        )
      },
      (None, None) => (msg_to_parent, PerEdge::new()),
      _ => {
        return Err(make_internal_report!(
          "Sparse marginal backward pass received an inconsistent parent edge at node {key}"
        ));
      },
    };

    let output = Arc::new((SparseNodeState { sequence, profile }, parent_message));
    Ok((Arc::clone(&output), output))
  }

  fn child_edges<'c>(
    &self,
    key: GraphNodeKey,
    children: &'c [SparseBackwardOutput],
  ) -> Result<Vec<&'c SparseEdgeBackward>, Report> {
    izip!(self.schedule.children(key), children)
      .map(|(&(child_key, edge_key), child)| {
        child
          .1
          .get(self.schedule.parent_slot(edge_key))
          .ok_or_else(|| make_internal_report!("Backward child {child_key} published no message for its parent edge"))
      })
      .collect()
  }

  fn combine_child_messages(
    &self,
    key: GraphNodeKey,
    child_edges: &[&SparseEdgeBackward],
    obs: &SparseNodeObs,
    sequence: &Seq,
  ) -> Result<SparseSeqDistribution, Report> {
    let partition = self.partition;
    let alphabet = &partition.alphabet;
    let child_keys = self.schedule.children(key);
    let child_messages = child_edges.iter().map(|edge| &edge.msg_from_child).collect_vec();
    let child_missing = child_edges
      .iter()
      .map(|edge| range_difference(&edge.child_no_data, &obs.non_char))
      .collect_vec();
    let n_events: usize = izip!(child_keys, &child_messages, &child_missing)
      .map(|(&(_, edge_key), message, missing)| {
        2 * partition.obs_edges[&edge_key].fitch_subs().len()
          + message.variable.len()
          + missing.iter().map(|(start, end)| end - start).sum::<usize>()
      })
      .sum();

    let mut events = CandidateSiteEvents::with_capacity(child_keys.len(), n_events);
    for (ci, (&(_, edge_key), child_message, missing)) in izip!(child_keys, &child_messages, &child_missing).enumerate()
    {
      for mutation in partition.obs_edges[&edge_key].fitch_subs() {
        events.insert(mutation.pos(), SiteSlot::State, mutation.reff());
        events.insert(mutation.pos(), SiteSlot::Reference(ci), mutation.qry());
      }
      for (pos, profile) in &child_message.variable {
        events.or_insert(*pos, SiteSlot::State, profile.state);
      }
      for pos in missing.iter().flat_map(|&(start, end)| start..end) {
        let state = sequence[pos];
        if alphabet.is_canonical(state) {
          events.or_insert(pos, SiteSlot::State, state);
        }
      }
    }
    let mut sites = events.finish();
    sites.fill_missing_references(|ci, pos, parent_state| {
      let child_obs = &partition.obs_nodes[&child_keys[ci].0];
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
      self.schedule.parents(key).is_empty().then_some(self.gtr.pi()),
    )
  }
}

type SparseBackwardOutput = Arc<(SparseNodeState, PerEdge<SparseEdgeBackward>)>;

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
