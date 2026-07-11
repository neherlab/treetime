use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::partition::sparse::{SparseEdgePartition, SparseNodePartition};
use crate::partition::traits::PartitionCompressed;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};
use treetime_primitives::seq;

#[derive(Clone, Debug)]
pub struct PartitionFitch {
  pub index: usize,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, SparseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, SparseEdgePartition>,
}

impl PartitionFitch {
  pub fn into_marginal_sparse<N, E>(self, gtr: GTR, graph: &Graph<N, E, ()>) -> Result<PartitionMarginalSparse, Report>
  where
    N: GraphNode,
    E: GraphEdge,
  {
    let root_key = graph.get_exactly_one_root()?.read_arc().key();
    let root_sequence = self.nodes[&root_key].seq.sequence.clone();
    let mut nodes = self.nodes;
    for (key, node_data) in &mut nodes {
      if *key != root_key && !graph.is_leaf(*key) {
        node_data.seq.sequence = seq![];
      }
    }
    Ok(PartitionMarginalSparse {
      index: self.index,
      gtr,
      alphabet: self.alphabet,
      length: self.length,
      root_sequence,
      nodes,
      edges: self.edges,
    })
  }

  pub fn into_marginal_dense(self, gtr: GTR) -> PartitionMarginalDense {
    PartitionMarginalDense::new(self.index, gtr, self.alphabet, self.length)
  }
}

impl PartitionCompressed for PartitionFitch {
  fn index(&self) -> usize {
    self.index
  }

  fn alphabet(&self) -> &Alphabet {
    &self.alphabet
  }

  fn length(&self) -> usize {
    self.length
  }

  fn nodes(&self) -> &BTreeMap<GraphNodeKey, SparseNodePartition> {
    &self.nodes
  }

  fn edges(&self) -> &BTreeMap<GraphEdgeKey, SparseEdgePartition> {
    &self.edges
  }

  fn nodes_mut(&mut self) -> &mut BTreeMap<GraphNodeKey, SparseNodePartition> {
    &mut self.nodes
  }

  fn edges_mut(&mut self) -> &mut BTreeMap<GraphEdgeKey, SparseEdgePartition> {
    &mut self.edges
  }

  fn storage_mut(
    &mut self,
  ) -> (
    &mut BTreeMap<GraphNodeKey, SparseNodePartition>,
    &mut BTreeMap<GraphEdgeKey, SparseEdgePartition>,
  ) {
    (&mut self.nodes, &mut self.edges)
  }
}
