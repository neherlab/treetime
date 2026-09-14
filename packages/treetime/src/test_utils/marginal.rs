use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::fitch::create_fitch_partition;
use crate::ancestral::marginal::profile_branch_lengths;
use crate::gtr::gtr::GTR;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::shared::update::{MarginalPasses, MarginalUpdate};
use crate::seq::alignment::get_common_length;
use eyre::Report;
use std::sync::LazyLock;
use treetime_graph::graph::Graph;
use treetime_io::fasta::read_many_fasta_str;
use treetime_io::nwk::{NwkParse, nwk_read_str};

pub static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

pub fn run_dense_marginal_with_newick(newick: &str, aln_str: &str, gtr: GTR) -> Result<f64, Report> {
  let NwkParse {
    graph,
    names,
    branch_lengths,
    ..
  } = nwk_read_str(newick)?;
  let graph: Graph = graph;
  let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;
  let length = get_common_length(&aln)?;
  let partition = PartitionMarginalDense::new(0, gtr, alphabet, length);

  let node_states = partition.attach_sequences(&graph, &aln, &names)?;
  let MarginalUpdate { log_lh, .. } =
    partition.marginal_update(&graph, &profile_branch_lengths(&branch_lengths), node_states)?;
  Ok(log_lh.value())
}

pub fn run_sparse_marginal_with_newick(newick: &str, aln_str: &str, gtr: GTR) -> Result<f64, Report> {
  let NwkParse {
    graph,
    names,
    branch_lengths,
    ..
  } = nwk_read_str(newick)?;
  let graph: Graph = graph;
  let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;

  let fitch = create_fitch_partition(&graph, 0, alphabet, &aln, &names)?;
  let (partition, node_states) = fitch.into_marginal_sparse(gtr, &graph)?;

  let MarginalUpdate { log_lh, .. } =
    partition.marginal_update(&graph, &profile_branch_lengths(&branch_lengths), node_states)?;
  Ok(log_lh.value())
}
