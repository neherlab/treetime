use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::fitch::create_fitch_partition;
use crate::ancestral::marginal::{initialize_marginal, marginal_update, profile_branch_lengths};
use crate::gtr::gtr::GTR;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::payload::ancestral::GraphAncestral;
use crate::seq::alignment::get_common_length;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::{Arc, LazyLock};
use treetime_io::fasta::read_many_fasta_str;
use treetime_io::nwk::{NwkParse, nwk_read_str};

pub static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

pub fn run_dense_marginal_with_newick(newick: &str, aln_str: &str, gtr: GTR) -> Result<f64, Report> {
  let NwkParse { graph, names, .. } = nwk_read_str(newick)?;
  let graph: GraphAncestral = graph;
  let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;
  let length = get_common_length(&aln)?;
  let partition = PartitionMarginalDense::new(0, gtr, alphabet, length);
  let partitions = [Arc::new(RwLock::new(partition))];

  initialize_marginal(&graph, &profile_branch_lengths(&graph), &partitions, &aln, &names).map(|log_lh| log_lh.value())
}

pub fn run_sparse_marginal_with_newick(newick: &str, aln_str: &str, gtr: GTR) -> Result<f64, Report> {
  let NwkParse { graph, names, .. } = nwk_read_str(newick)?;
  let graph: GraphAncestral = graph;
  let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;

  let fitch = create_fitch_partition(&graph, 0, alphabet, &aln, &names)?;
  let partition = fitch.into_marginal_sparse(gtr, &graph)?;
  let partitions = [Arc::new(RwLock::new(partition))];

  marginal_update(&graph, &profile_branch_lengths(&graph), &partitions).map(|log_lh| log_lh.value())
}
