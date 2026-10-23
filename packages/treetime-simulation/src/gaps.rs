use crate::evolve::validate_edge_duration;
use eyre::Report;
use itertools::Itertools;
use rand::Rng;
use rand::distr::weighted::WeightedIndex;
use rand::seq::index::sample;
use rand_distr::Distribution;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::{make_error, make_internal_error, make_report};

pub const GAP: AsciiChar = AsciiChar::from_byte_unchecked(b'-');

pub fn disjoint_gap_events(
  durations: &BTreeMap<GraphEdgeKey, f64>,
  count: usize,
  gap_length: usize,
  sequence_length: usize,
  rng: &mut impl Rng,
) -> Result<Vec<GapEvent>, Report> {
  if gap_length == 0 {
    return make_error!("Gap length must be positive, got 0");
  }
  let slot_width = gap_length + 1;
  #[expect(
    clippy::integer_division,
    reason = "the slot count is the floor of the sites over the slot width"
  )]
  let slot_count = (sequence_length + 1) / slot_width;
  if count > slot_count {
    return make_error!(
      "A sequence of {sequence_length} sites holds at most {slot_count} gaps of {gap_length} sites with a site between them, got {count} gaps"
    );
  }
  for (&edge, &duration) in durations {
    validate_edge_duration(edge, duration)?;
  }
  let edges = durations.keys().copied().collect_vec();
  let edge_choice = WeightedIndex::new(durations.values().copied())
    .map_err(|error| make_report!("Gap events need an edge with a positive duration: {error}"))?;
  let starts = sample(rng, slot_count, count)
    .into_iter()
    .map(|slot| slot * slot_width)
    .sorted()
    .collect_vec();
  Ok(
    starts
      .into_iter()
      .map(|start| GapEvent {
        edge: edges[edge_choice.sample(rng)],
        start,
        length: gap_length,
      })
      .collect(),
  )
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GapEvent {
  pub edge: GraphEdgeKey,
  pub start: usize,
  pub length: usize,
}

impl GapEvent {
  pub fn end(&self) -> usize {
    self.start + self.length
  }
}

pub fn apply_gap_events(
  graph: &Graph,
  sequences: &BTreeMap<GraphNodeKey, Seq>,
  events: &[GapEvent],
) -> Result<BTreeMap<GraphNodeKey, Seq>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  let events_by_edge = events.iter().into_group_map_by(|event| event.edge);
  for event in events {
    let Some(edge) = graph.get_edge(event.edge) else {
      return make_error!(
        "Gap event at site {} is on edge {}, which the tree has not",
        event.start,
        event.edge
      );
    };
    let length = sequences
      .get(&edge.target())
      .ok_or_else(|| make_report!("Node {} has no sequence", edge.target()))?
      .len();
    if event.length == 0 || event.end() > length {
      return make_error!(
        "Gap event on edge {} covers sites {}..{}, outside the {length} sites of the sequence",
        event.edge,
        event.start,
        event.end()
      );
    }
  }

  let mut gapped: BTreeMap<GraphNodeKey, Seq> = BTreeMap::new();
  for &node in schedule.preorder() {
    let Some(seq) = sequences.get(&node) else {
      return make_error!("Node {node} has no sequence");
    };
    let mut seq = seq.clone();
    match schedule.parents(node) {
      [] => {},
      [(parent, edge)] => {
        let parent_seq = &gapped[parent];
        if parent_seq.len() != seq.len() {
          return make_error!(
            "Sequences of nodes {parent} and {node} have different lengths: {} and {}",
            parent_seq.len(),
            seq.len()
          );
        }
        for (site, &parent_char) in seq.as_mut_slice().iter_mut().zip(parent_seq.iter()) {
          if parent_char == GAP {
            *site = GAP;
          }
        }
        for event in events_by_edge.get(edge).into_iter().flatten() {
          seq.as_mut_slice()[event.start..event.end()].fill(GAP);
        }
      },
      parents => {
        return make_internal_error!("Node {node} has {} parents in a tree", parents.len());
      },
    }
    gapped.insert(node, seq);
  }
  Ok(gapped)
}
