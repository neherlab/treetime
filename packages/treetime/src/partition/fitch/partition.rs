use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::sparse::{
  FitchNodeData, SparseEdgeObs, SparseNodeObs, SparseNodeState, SparseSeqDistribution,
};
use crate::seq::mutation::Sub;
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AsciiChar, Seq, seq};

#[derive(Clone, Debug, Serialize)]
pub struct PartitionFitch {
  pub index: usize,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, FitchNodeData>,
  pub edges: BTreeMap<GraphEdgeKey, SparseEdgeObs>,
}

impl PartitionFitch {
  /// Hand the Fitch results off to the sparse marginal representation: split each node's Fitch data
  /// into the durable observations the partition owns and the seed node state the marginal passes
  /// evolve. Internal-node sequences are cleared (the forward pass rebuilds them); leaf and root
  /// sequences are kept.
  pub fn into_marginal_sparse(
    self,
    graph: &Graph,
  ) -> Result<(PartitionMarginalSparse, BTreeMap<GraphNodeKey, SparseNodeState>), Report> {
    let root_key = graph.get_exactly_one_root()?.read_arc().key();
    let root_sequence = self.nodes[&root_key].seq.sequence.clone();

    let mut obs_nodes = BTreeMap::new();
    let mut node_states = BTreeMap::new();
    for (key, node_data) in self.nodes {
      let keep_sequence = key == root_key || graph.is_leaf(key);
      let sequence = if keep_sequence { node_data.seq.sequence } else { seq![] };
      node_states.insert(
        key,
        SparseNodeState {
          sequence,
          profile: SparseSeqDistribution::default(),
          emitted: None,
        },
      );
      obs_nodes.insert(
        key,
        SparseNodeObs {
          unknown: node_data.seq.unknown,
          gaps: node_data.seq.gaps,
          non_char: node_data.seq.non_char,
          composition: node_data.seq.composition,
          fitch: node_data.seq.fitch,
        },
      );
    }

    let partition = PartitionMarginalSparse {
      index: self.index,
      alphabet: self.alphabet,
      length: self.length,
      root_sequence,
      obs_nodes,
      obs_edges: self.edges,
    };
    Ok((partition, node_states))
  }

  pub fn into_marginal_dense(self) -> PartitionMarginalDense {
    PartitionMarginalDense::new(self.index, self.alphabet, self.length)
  }

  pub fn sequence_length(&self) -> usize {
    self.length
  }

  pub fn edge_subs(&self, _graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    Ok(self.edges[&edge_key].fitch_subs().to_vec())
  }

  pub fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<crate::seq::indel::InDel> {
    self.edges[&edge_key].indels.clone()
  }

  pub fn root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    Ok(self.nodes[&graph.root_key()?].seq.sequence.clone())
  }

  pub fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    self.nodes[&node_key].seq.sequence.clone()
  }

  pub fn edge_effective_length(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
    Ok(
      self.nodes[&parent_key]
        .seq
        .sequence
        .iter()
        .zip(&self.nodes[&child_key].seq.sequence)
        .filter(|(parent, child)| self.alphabet.is_canonical(**parent) && self.alphabet.is_canonical(**child))
        .count(),
    )
  }

  /// The alphabet's ambiguous (unknown) character, used to fill masked positions in output sequences.
  pub fn ambiguous_char(&self) -> AsciiChar {
    self.alphabet.unknown()
  }
}
