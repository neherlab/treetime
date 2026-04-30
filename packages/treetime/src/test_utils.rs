use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::commands::ancestral::fitch::get_common_length;
use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::gtr::gtr::GTR;
use crate::representation::partition::fitch::PartitionFitch;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::{Arc, LazyLock};
use treetime_graph::edge::GraphEdge;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node::Named;
use treetime_io::fasta::read_many_fasta_str;
use treetime_io::nwk::nwk_read_str;

/// Default nucleotide alphabet for test helpers.
pub static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

/// Run dense marginal reconstruction on a tree given as Newick string.
///
/// Returns the log-likelihood of the reconstruction.
pub fn run_dense_marginal_with_newick(newick: &str, aln_str: &str, gtr: GTR) -> Result<f64, Report> {
  let graph: GraphAncestral = nwk_read_str(newick)?;
  let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;
  let length = get_common_length(&aln)?;
  let partition = PartitionMarginalDense::new(0, gtr, alphabet, length);
  let partitions = [Arc::new(RwLock::new(partition))];

  initialize_marginal(&graph, &partitions, &aln)
}

/// Run sparse marginal reconstruction on a tree given as Newick string.
///
/// Returns the log-likelihood of the reconstruction.
pub fn run_sparse_marginal_with_newick(newick: &str, aln_str: &str, gtr: GTR) -> Result<f64, Report> {
  let graph: GraphAncestral = nwk_read_str(newick)?;
  let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;

  let fitch = PartitionFitch::compress(&graph, 0, alphabet, &aln)?;
  let partition = fitch.into_marginal_sparse(gtr, &graph)?;
  let partitions = [Arc::new(RwLock::new(partition))];

  update_marginal(&graph, &partitions)
}

/// Find a node key by its name in a graph where nodes implement `Named`.
pub fn find_node_key_by_name<N, E, D>(graph: &Graph<N, E, D>, name: &str) -> Option<GraphNodeKey>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Send + Sync,
{
  for node in graph.get_nodes() {
    let node = node.read_arc();
    let payload = node.payload().read_arc();
    if payload.name().is_some_and(|n| n.as_ref() == name) {
      return Some(node.key());
    }
  }
  None
}

/// Find an edge key by source and target node names.
pub fn find_edge_key<N, E, D>(graph: &Graph<N, E, D>, source_name: &str, target_name: &str) -> Option<GraphEdgeKey>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Send + Sync,
{
  let source_key = find_node_key_by_name(graph, source_name)?;
  let target_key = find_node_key_by_name(graph, target_name)?;

  for edge in graph.get_edges() {
    let edge = edge.read_arc();
    if edge.source() == source_key && edge.target() == target_key {
      return Some(edge.key());
    }
  }
  None
}
