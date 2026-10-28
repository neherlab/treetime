use crate::partition::marginal::reconstruction::SparseReconstruction;
use crate::partition::marginal::sparse::reconstruct::parsimony_seq;
use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs};
use crate::seq::composition::Composition;
use crate::seq::indel::compute_node_ranges;
use eyre::Report;
use itertools::Itertools;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;
use treetime_utils::interval::range_difference::range_difference;
use treetime_utils::make_internal_report;

pub(crate) fn new_node_obs(
  graph: &Graph,
  reconstruction: &SparseReconstruction,
  edge_key: GraphEdgeKey,
  child_keys: &[GraphNodeKey],
) -> Result<SparseNodeObs, Report> {
  let partition = &reconstruction.partition;
  let alphabet = &partition.alphabet;
  let children = child_keys.iter().map(|key| &partition.obs_nodes[key]).collect_vec();
  let child_non_chars = children.iter().map(|child| &child.non_char).collect_vec();
  let child_gaps = children.iter().map(|child| &child.gaps).collect_vec();
  let ranges = compute_node_ranges(&child_non_chars, &child_gaps);
  let gaps = range_difference(&ranges.non_char, &ranges.unknown);

  let mut sequence = node_sequence(graph, reconstruction, graph.get_source_node_key(edge_key)?)?;
  for sub in partition.obs_edges[&edge_key].fitch_subs() {
    sequence[sub.pos()] = sub.qry();
  }
  for &(start, end) in &gaps {
    sequence[start..end].fill(alphabet.gap());
  }
  for &(start, end) in &ranges.unknown {
    sequence[start..end].fill(alphabet.unknown());
  }

  let mut obs = SparseNodeObs::empty(alphabet);
  obs.composition = Composition::with_seq(&sequence, alphabet);
  obs.unknown = ranges.unknown;
  obs.gaps = gaps;
  obs.non_char = ranges.non_char;
  Ok(obs)
}

fn node_sequence(graph: &Graph, reconstruction: &SparseReconstruction, node_key: GraphNodeKey) -> Result<Seq, Report> {
  let partition = &reconstruction.partition;
  let mut path = Vec::new();
  let mut current = node_key;
  let base = loop {
    if let Some(state) = reconstruction.node_states.get(current)
      && !state.sequence.is_empty()
    {
      break (*state.sequence).clone();
    }
    let node = graph
      .get_node(current)
      .ok_or_else(|| make_internal_report!("Node {current} not found"))?;
    match node.inbound() {
      [] => break partition.root_sequence.clone(),
      &[edge_key] => {
        path.push((edge_key, current));
        current = graph.get_source_node_key(edge_key)?;
      },
      _ => return Err(make_internal_report!("Node {current} has more than one parent")),
    }
  };
  let no_mutations = SparseEdgeObs::default();
  Ok(path.iter().rev().fold(base, |sequence, &(edge_key, key)| {
    parsimony_seq(
      &sequence,
      partition.obs_edges.get(edge_key).unwrap_or(&no_mutations),
      &partition.obs_nodes[&key],
      &partition.alphabet,
    )
  }))
}
