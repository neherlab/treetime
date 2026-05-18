use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::fitch::create_fitch_partition;
use crate::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::gtr::gtr::GTR;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::payload::ancestral::GraphAncestral;
use crate::seq::alignment::get_common_length;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::{Arc, LazyLock};
use treetime_io::fasta::read_many_fasta_str;
use treetime_io::nwk::nwk_read_str;

pub static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

pub fn run_dense_marginal_with_newick(newick: &str, aln_str: &str, gtr: GTR) -> Result<f64, Report> {
  let graph: GraphAncestral = nwk_read_str(newick)?;
  let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;
  let length = get_common_length(&aln)?;
  let partition = PartitionMarginalDense::new(0, gtr, alphabet, length);
  let partitions = [Arc::new(RwLock::new(partition))];

  initialize_marginal(&graph, &partitions, &aln)
}

pub fn run_sparse_marginal_with_newick(newick: &str, aln_str: &str, gtr: GTR) -> Result<f64, Report> {
  let graph: GraphAncestral = nwk_read_str(newick)?;
  let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;

  let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;
  let partition = fitch.into_marginal_sparse(gtr, &graph)?;
  let partitions = [Arc::new(RwLock::new(partition))];

  update_marginal(&graph, &partitions)
}
