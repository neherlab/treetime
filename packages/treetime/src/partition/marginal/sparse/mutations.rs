use crate::alphabet::alphabet::Alphabet;
use crate::partition::marginal::sample::{Resolve, resolve_profile};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::marginal::sparse::reconstruct::impute_state;
use crate::partition::storage::sparse::{SparseEdgeForward, SparseNodeObs, SparseNodeState, SparseSeqDistribution};
use crate::seq::mutation::{Mutation, MutationTrack, Sub, combine_edge_mutations};
use crate::{make_internal_error, make_internal_report};
use eyre::Report;
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::AsciiChar;
use treetime_utils::interval::range_union::range_union;

pub(crate) fn sparse_edge_mutations(
  partition: &PartitionMarginalSparse,
  graph: &Graph,
  node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
  forward: &BTreeMap<GraphEdgeKey, SparseEdgeForward>,
  impute: bool,
  track: &MutationTrack,
) -> Result<BTreeMap<GraphEdgeKey, Vec<Mutation>>, Report> {
  let edges = graph
    .get_nodes()
    .map(|node| Ok(graph.node_parent(node.key())?.map(|parent| (node.key(), parent))))
    .filter_map_ok(|edge| edge)
    .collect::<Result<Vec<_>, Report>>()?;
  edges
    .into_par_iter()
    .map(|(child_key, (parent_key, edge_key))| {
      let child = EdgeChild::new(
        partition,
        graph,
        node_states,
        forward,
        impute,
        child_key,
        parent_key,
        edge_key,
      )?;
      let subs = edge_sequence_subs(partition, &node_states[&parent_key], &child, edge_key)?;
      let mutations = combine_edge_mutations(subs, &partition.obs_edges[&edge_key].indels, track)?;
      Ok((edge_key, mutations))
    })
    .collect()
}

fn edge_sequence_subs(
  partition: &PartitionMarginalSparse,
  parent: &SparseNodeState,
  child: &EdgeChild<'_>,
  edge_key: GraphEdgeKey,
) -> Result<Vec<Sub>, Report> {
  let alphabet = &partition.alphabet;
  if parent.sequence.len() != child.state.sequence.len() {
    return make_internal_error!(
      "Parent sequence has length {}, but child sequence has length {}",
      parent.sequence.len(),
      child.state.sequence.len()
    );
  }
  let edge_obs = &partition.obs_edges[&edge_key];

  let mut positions: BTreeSet<usize> = parent.profile.variable.keys().copied().collect();
  positions.extend(edge_obs.fitch_subs().iter().map(Sub::pos));
  if child.leaf.is_some() {
    positions.extend(child.obs.fitch.variable.keys());
  } else {
    positions.extend(child.state.profile.variable.keys());
  }
  let ranges = range_union(&[
    edge_obs.indels.iter().map(|indel| indel.range).collect_vec(),
    child.obs.unknown.clone(),
  ]);

  positions
    .iter()
    .copied()
    .merge(ranges.iter().flat_map(|&(start, end)| start..end))
    .dedup()
    .filter_map(|pos| {
      let reff = internal_state(parent, pos, alphabet);
      let qry = child.state_at(pos, alphabet);
      (reff != qry && !alphabet.is_gap(reff) && !alphabet.is_gap(qry)).then(|| Sub::new(reff, pos, qry))
    })
    .collect()
}

fn internal_state(node: &SparseNodeState, pos: usize, alphabet: &Alphabet) -> AsciiChar {
  let state = node.sequence[pos];
  if state == alphabet.gap() {
    return state;
  }
  node.profile.variable.get(&pos).map_or(state, |var| {
    alphabet.char(resolve_profile(var.dis.view(), &mut Resolve::Argmax))
  })
}

struct EdgeChild<'a> {
  state: &'a SparseNodeState,
  obs: &'a SparseNodeObs,
  leaf: Option<LeafImputation<'a>>,
}

impl<'a> EdgeChild<'a> {
  #[expect(
    clippy::too_many_arguments,
    reason = "the child of an edge is located by its graph, its partition and its pass results"
  )]
  fn new(
    partition: &'a PartitionMarginalSparse,
    graph: &Graph,
    node_states: &'a BTreeMap<GraphNodeKey, SparseNodeState>,
    forward: &'a BTreeMap<GraphEdgeKey, SparseEdgeForward>,
    impute: bool,
    child_key: GraphNodeKey,
    parent_key: GraphNodeKey,
    edge_key: GraphEdgeKey,
  ) -> Result<Self, Report> {
    let is_leaf = graph
      .get_node(child_key)
      .ok_or_else(|| make_internal_report!("Node {child_key} not found while deriving its mutations"))?
      .is_leaf();
    let leaf = is_leaf.then(|| LeafImputation {
      impute: impute.then(|| (&forward[&edge_key].msg_from_parent, &node_states[&parent_key])),
    });
    Ok(Self {
      state: &node_states[&child_key],
      obs: &partition.obs_nodes[&child_key],
      leaf,
    })
  }

  fn state_at(&self, pos: usize, alphabet: &Alphabet) -> AsciiChar {
    match &self.leaf {
      None => internal_state(self.state, pos, alphabet),
      Some(leaf) => leaf.state_at(self.state, self.obs, pos, alphabet),
    }
  }
}

struct LeafImputation<'a> {
  impute: Option<(&'a SparseSeqDistribution, &'a SparseNodeState)>,
}

impl LeafImputation<'_> {
  fn state_at(&self, node: &SparseNodeState, node_obs: &SparseNodeObs, pos: usize, alphabet: &Alphabet) -> AsciiChar {
    let ambiguous = node_obs.fitch.variable.get(&pos);
    let state = ambiguous.map_or(node.sequence[pos], |states| alphabet.set_to_char(*states));
    let Some((down, parent)) = self.impute else {
      return state;
    };
    let visits = usize::from(sorted_ranges_contain(&node_obs.unknown, pos)) + usize::from(ambiguous.is_some());
    (0..visits).fold(state, |state, _| impute_state(state, pos, down, parent, alphabet))
  }
}

fn sorted_ranges_contain(ranges: &[(usize, usize)], pos: usize) -> bool {
  let index = ranges.partition_point(|&(_, end)| end <= pos);
  ranges.get(index).is_some_and(|&(start, _)| start <= pos)
}
