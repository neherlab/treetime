use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::fitch::create_fitch_partition;
use crate::ancestral::marginal::branch_lengths_or_zero;
use crate::gtr::gtr::GTR;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::shared::update::{MarginalPasses, MarginalUpdate};
use crate::seq::alignment::get_common_length;
use eyre::Report;
use std::sync::LazyLock;
use treetime_graph::graph::Graph;
use treetime_io::fasta::read_many_fasta_str;
use treetime_io::nwk::nwk_fasta_node_inputs;
use treetime_io::nwk::nwk_read_str;

pub static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

pub fn run_dense_marginal_with_newick(newick: &str, aln_str: &str, gtr: &GTR) -> Result<f64, Report> {
  let nwk_parsed = nwk_read_str(newick)?;
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let graph: Graph = graph;
  let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;
  let length = get_common_length(&aln)?;
  let partition = PartitionMarginalDense::new(0, alphabet, length);

  let node_states = partition.attach_sequences(&graph, &nwk_fasta_node_inputs(&graph, &names, aln))?;
  let MarginalUpdate { log_lh, .. } =
    partition.marginal_update(gtr, &graph, &branch_lengths_or_zero(&branch_lengths), node_states)?;
  Ok(log_lh.value())
}

pub fn run_sparse_marginal_with_newick(newick: &str, aln_str: &str, gtr: &GTR) -> Result<f64, Report> {
  let nwk_parsed = nwk_read_str(newick)?;
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let graph: Graph = graph;
  let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;

  let fitch = create_fitch_partition(&graph, 0, alphabet, &nwk_fasta_node_inputs(&graph, &names, aln))?;
  let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;

  let MarginalUpdate { log_lh, .. } =
    partition.marginal_update(gtr, &graph, &branch_lengths_or_zero(&branch_lengths), node_states)?;
  Ok(log_lh.value())
}
