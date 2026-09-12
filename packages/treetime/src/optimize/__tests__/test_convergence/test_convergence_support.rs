#[cfg(test)]
pub mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::{initialize_marginal, marginal_update, profile_branch_lengths};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::optimize::run_loop::optimize_partition_view;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use indoc::indoc;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;

  use std::sync::LazyLock;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};

  pub static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  // Small tree with 4 leaves
  pub const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  pub fn simple_alignment() -> Result<Vec<FastaRecord>, Report> {
    read_many_fasta_str(
      indoc! {r#"
      >A
      ACGTACGTACGTACGT
      >B
      ACGTACGTACGTACGA
      >C
      ACGTACGTACGTACGG
      >D
      ACGTACGTACGTACGC
    "#},
      &*NUC_ALPHABET,
    )
  }

  pub fn setup_partitions(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[FastaRecord],
    branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<(Vec<PartitionMarginalDense>, Vec<PartitionMarginalSparse>), Report> {
    let alphabet_dense = Alphabet::new(AlphabetName::Nuc)?;
    let alphabet_sparse = Alphabet::new(AlphabetName::Nuc)?;

    let mut dense_partitions = vec![PartitionMarginalDense::new(
      0,
      jc69(JC69Params::default())?,
      alphabet_dense,
      get_common_length(aln)?,
    )];

    let fitch = create_fitch_partition(graph, 1, alphabet_sparse, aln, names)?;
    let mut sparse_partitions = vec![fitch.into_marginal_sparse(jc69(JC69Params::default())?, graph)?];
    initialize_marginal(
      graph,
      &profile_branch_lengths(branch_lengths),
      &mut dense_partitions,
      aln,
      names,
    )?
    .value();
    marginal_update(graph, &profile_branch_lengths(branch_lengths), &mut sparse_partitions)?.value();

    let mixed_partitions = optimize_partition_view(&dense_partitions, &sparse_partitions);
    initial_guess_mixed(graph, &mixed_partitions, true, false, branch_lengths)?;

    Ok((dense_partitions, sparse_partitions))
  }

  pub fn compute_total_lh(
    graph: &Graph,
    dense_partitions: &mut [PartitionMarginalDense],
    sparse_partitions: &mut [PartitionMarginalSparse],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<f64, Report> {
    let dense_lh = marginal_update(graph, &profile_branch_lengths(branch_lengths), dense_partitions)?.value();
    let sparse_lh = marginal_update(graph, &profile_branch_lengths(branch_lengths), sparse_partitions)?.value();
    Ok(dense_lh + sparse_lh)
  }
}
