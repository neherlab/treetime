use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::partition::marginal::shared::update::MarginalBackward;
use crate::partition::marginal::sparse::candidate_sites::{CandidateSiteEvents, SiteSlot};
use crate::partition::marginal::sparse::globally_missing::GloballyMissing;
use crate::partition::marginal::sparse::message::{ScaleInputs, combine_messages, propagate_raw};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::marginal::sparse::sequences::fitch_node_sequences;
use crate::partition::storage::sparse::{
  FixedProfiles, SiteScales, SparseEdgeBackward, SparseEdgeObs, SparseNodeObs, SparseNodeState, SparseSeqDistribution,
};
use crate::partition::storage::var_pos_map::VarPosMap;
use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::aview_mut1;
use smallvec::smallvec;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::schedule_scan::{NodeEdgeMaps, PerEdge};
use treetime_primitives::{AsciiChar, LogLh, Seq};
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
  let inheritance = &partition.inheritance;
  inheritance.require_one_parent(&schedule)?;
  let sequences = fitch_node_sequences(partition, &schedule, node_states)?;
  let pass = SparseBackwardPass {
    schedule: &schedule,
    partition,
    gtr,
    branch_lengths,
  };
  let edges = partition.obs_edges.iter().collect();
  let outputs = schedule.scan_up(
    inheritance.join_parent_edge(&schedule, sequences, edges),
    |key, (sequence, edge_obs), children| match sequence {
      Some(sequence) => pass.node(key, sequence, edge_obs, &children),
      None => Ok((None, None)),
    },
  )?;
  let NodeEdgeMaps { nodes, edges } = inheritance.split_parent_edge(&schedule, outputs.unwrap_present().into_map()?)?;
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

    let inheritance = &partition.inheritance;
    let children = inheritance.child_messages(self.schedule, key, children)?;
    let (msg_to_parent, no_data) = if children.is_empty() {
      (
        leaf_message(obs, &partition.alphabet)?,
        range_union(&[obs.gaps.clone(), obs.unknown.clone()]),
      )
    } else {
      let child_edges = child_edges(key, &children)?;
      let no_data = range_intersection_iter(child_edges.iter().map(|edge| &edge.child_no_data)).collect_vec();
      let msg = self.combine_child_messages(key, &children, &child_edges, obs, &sequence, &no_data)?;
      (msg, no_data)
    };

    let (profile, globally_missing, parent_message) = match (inheritance.parent(self.schedule, key), edge_obs) {
      (Some((_, _, edge_key)), Some(_)) => {
        let branch_length = fix_branch_length(partition.total_length, self.branch_lengths[&edge_key]);
        let msg_from_child = propagate_raw(&gtr.expQt(branch_length).t().to_owned(), &msg_to_parent);
        (
          SparseSeqDistribution::default(),
          Arc::default(),
          smallvec![SparseEdgeBackward {
            msg_to_parent,
            msg_from_child,
            child_no_data: no_data,
          }],
        )
      },
      (None, None) => (
        msg_to_parent,
        Arc::new(GloballyMissing::new(no_data, gtr.pi(), &partition.alphabet)),
        PerEdge::new(),
      ),
      _ => {
        return Err(make_internal_report!(
          "Sparse marginal backward pass received an inconsistent parent edge at node {key}"
        ));
      },
    };

    let output = Arc::new((
      SparseNodeState {
        sequence,
        profile,
        globally_missing,
      },
      parent_message,
    ));
    Ok((Some(Arc::clone(&output)), Some(output)))
  }

  fn combine_child_messages(
    &self,
    key: GraphNodeKey,
    children: &[CarriedChild<'_>],
    child_edges: &[&SparseEdgeBackward],
    obs: &SparseNodeObs,
    sequence: &Seq,
    no_data: &[(usize, usize)],
  ) -> Result<SparseSeqDistribution, Report> {
    let partition = self.partition;
    let alphabet = &partition.alphabet;
    let child_keys = children
      .iter()
      .map(|&(child_key, edge_key, _)| (child_key, edge_key))
      .collect_vec();
    let child_messages = child_edges.iter().map(|edge| &edge.msg_from_child).collect_vec();
    let child_missing = child_edges
      .iter()
      .map(|edge| range_difference(&edge.child_no_data, &obs.non_char))
      .collect_vec();
    let n_events: usize = izip!(&child_keys, &child_messages, &child_missing)
      .map(|(&(_, edge_key), message, missing)| {
        2 * partition.obs_edges[&edge_key].fitch_subs().len()
          + message.variable.len()
          + missing.iter().map(|(start, end)| end - start).sum::<usize>()
      })
      .sum();

    let gapped = range_difference(&obs.non_char, no_data);
    let n_events = n_events + gapped.iter().map(|(start, end)| end - start).sum::<usize>();
    let mut events = CandidateSiteEvents::with_capacity(child_keys.len(), n_events);
    for pos in gapped.iter().flat_map(|&(start, end)| start..end) {
      events.or_insert(pos, SiteSlot::State, alphabet.gap());
    }
    for (ci, (&(_, edge_key), child_message, missing)) in
      izip!(&child_keys, &child_messages, &child_missing).enumerate()
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
    sites.relabel(|pos, state| own_state(alphabet, obs, sequence, pos, state));
    sites.fill_missing_references(|ci, pos, parent_state| {
      if range_contains(&child_edges[ci].child_no_data, pos) {
        if range_contains(&partition.obs_nodes[&child_keys[ci].0].gaps, pos) {
          alphabet.gap()
        } else {
          alphabet.unknown()
        }
      } else {
        let child_state = children[ci].2.0.sequence[pos];
        if alphabet.is_canonical(child_state) {
          child_state
        } else {
          parent_state
        }
      }
    });

    let child_sequences = children.iter().map(|child| &*child.2.0.sequence).collect_vec();
    combine_messages(
      &obs.composition,
      &child_messages,
      &sites,
      alphabet,
      partition
        .inheritance
        .parent(self.schedule, key)
        .is_none()
        .then_some(self.gtr.pi()),
      partition.track_site_scales.then_some(ScaleInputs {
        sequence,
        child_sequences: &child_sequences,
      }),
    )
  }
}

type SparseBackwardOutput = Option<Arc<(SparseNodeState, PerEdge<SparseEdgeBackward>)>>;

type CarriedChild<'c> = (
  GraphNodeKey,
  GraphEdgeKey,
  &'c Arc<(SparseNodeState, PerEdge<SparseEdgeBackward>)>,
);

fn child_edges<'c>(key: GraphNodeKey, children: &[CarriedChild<'c>]) -> Result<Vec<&'c SparseEdgeBackward>, Report> {
  children
    .iter()
    .map(|&(child_key, _, child)| {
      child.1.first().ok_or_else(|| {
        make_internal_report!("Backward child {child_key} of node {key} published no message for its parent edge")
      })
    })
    .collect()
}

pub(crate) fn own_state(
  alphabet: &Alphabet,
  obs: &SparseNodeObs,
  sequence: &Seq,
  pos: usize,
  state: AsciiChar,
) -> AsciiChar {
  if range_contains(&obs.non_char, pos) {
    alphabet.gap()
  } else if alphabet.is_canonical(sequence[pos]) {
    sequence[pos]
  } else {
    state
  }
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
    scales: SiteScales::default(),
  })
}
