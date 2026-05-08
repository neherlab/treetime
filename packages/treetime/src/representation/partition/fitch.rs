use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch::compress_sequences;
use crate::gtr::get_gtr::{GtrModelName, get_gtr_by_name};
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::fitch::infer_gtr_fitch;
use crate::make_report;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition::traits::PartitionCompressed;
use crate::representation::payload::sparse::{SparseEdgePartition, SparseNodePartition};
use crate::seq::alignment::get_common_length;
use eyre::Report;
use maplit::btreemap;
use parking_lot::RwLock;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey, NodeAncestralOps};
use treetime_io::fasta::FastaRecord;
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
  /// Phase 1: construct and compress in one step.
  pub fn compress<N, E>(
    graph: &Graph<N, E, ()>,
    index: usize,
    alphabet: Alphabet,
    aln: &[FastaRecord],
  ) -> Result<Self, Report>
  where
    N: NodeAncestralOps,
    E: GraphEdge,
  {
    let length = get_common_length(aln)?;
    let partition = Arc::new(RwLock::new(Self {
      index,
      alphabet,
      length,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));
    compress_sequences(graph, std::slice::from_ref(&partition), aln)?;
    Arc::try_unwrap(partition)
      .map(|rw| rw.into_inner())
      .map_err(|_arc| make_report!("PartitionFitch::compress: Arc still shared after compress_sequences"))
  }

  /// Phase 2: infer GTR from Fitch substitution counts.
  pub fn infer_gtr<N, E, D>(&self, graph: &Graph<N, E, D>) -> Result<GTR, Report>
  where
    N: GraphNode,
    E: GraphEdge + HasBranchLength,
    D: Send + Sync,
  {
    infer_gtr_fitch(self, graph)
  }

  /// Phase 2 (dispatch): infer or resolve named model.
  pub fn resolve_gtr<N, E, D>(&self, graph: &Graph<N, E, D>, model_name: GtrModelName) -> Result<GTR, Report>
  where
    N: GraphNode,
    E: GraphEdge + HasBranchLength,
    D: Send + Sync,
  {
    match model_name {
      GtrModelName::Infer => self.infer_gtr(graph),
      _ => get_gtr_by_name(model_name),
    }
  }

  /// Phase 3 (sparse): move Fitch data into marginal partition (zero-copy for fields).
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

  /// Phase 3 (dense): take config, discard Fitch data.
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
}
