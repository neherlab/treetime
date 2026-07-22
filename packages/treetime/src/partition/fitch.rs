use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::partition::sparse::{SparseEdgePartition, SparseNodePartition};
use crate::partition::traits::{BranchTopology, PartitionBranchOps, PartitionCompressed};
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};
use treetime_primitives::{Seq, seq};

#[derive(Clone, Debug, Serialize)]
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

impl PartitionBranchOps for PartitionFitch {
  fn sequence_length(&self) -> usize {
    self.length
  }

  fn edge_subs(
    &self,
    _graph: &dyn BranchTopology,
    edge_key: GraphEdgeKey,
  ) -> Result<Vec<crate::seq::mutation::Sub>, Report> {
    Ok(self.edges[&edge_key].fitch_subs().to_vec())
  }

  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<crate::seq::indel::InDel> {
    self.edges[&edge_key].indels.clone()
  }

  fn root_sequence(&self, graph: &dyn BranchTopology) -> Result<Seq, Report> {
    Ok(self.nodes[&graph.root_key()?].seq.sequence.clone())
  }

  fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    self.nodes[&node_key].seq.sequence.clone()
  }

  fn edge_effective_length(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<usize, Report> {
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
}
