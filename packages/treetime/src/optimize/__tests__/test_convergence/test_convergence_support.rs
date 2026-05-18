#[cfg(test)]
pub mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::optimize::run_loop::collect_optimize_partitions;
  use crate::partition::fitch::PartitionFitch;
  use crate::partition::marginal_dense::PartitionMarginalDense;
  use crate::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::partition::traits::PartitionOptimizeVec;
  use crate::partition::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::sync::{Arc, LazyLock};
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
    graph: &GraphAncestral,
    aln: &[FastaRecord],
  ) -> Result<
    (
      Vec<Arc<RwLock<PartitionMarginalDense>>>,
      Vec<Arc<RwLock<PartitionMarginalSparse>>>,
      PartitionOptimizeVec,
    ),
    Report,
  > {
    let alphabet_dense = Alphabet::new(AlphabetName::Nuc)?;
    let alphabet_sparse = Alphabet::new(AlphabetName::Nuc)?;

    let dense_partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: alphabet_dense,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    let fitch = PartitionFitch::compress(graph, 1, alphabet_sparse, aln)?;
    let sparse_partitions = vec![Arc::new(RwLock::new(
      fitch.into_marginal_sparse(jc69(JC69Params::default())?, graph)?,
    ))];
    initialize_marginal(graph, &dense_partitions, aln)?;
    update_marginal(graph, &sparse_partitions)?;

    let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);
    initial_guess_mixed(graph, &mixed_partitions, true)?;

    Ok((dense_partitions, sparse_partitions, mixed_partitions))
  }

  pub fn compute_total_lh(
    graph: &GraphAncestral,
    dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
    sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  ) -> Result<f64, Report> {
    let dense_lh = update_marginal(graph, dense_partitions)?;
    let sparse_lh = update_marginal(graph, sparse_partitions)?;
    Ok(dense_lh + sparse_lh)
  }
}
