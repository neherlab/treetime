#[cfg(test)]
pub mod tests {
  //! Tests for dense vs sparse optimization equivalence.
  //!
  //! These tests verify that dense and sparse optimization modes produce consistent
  //! results. Due to implementation differences (different zero-branch thresholds,
  //! optimization paths), the modes may not produce identical results but should:
  //! 1. Both produce finite, valid log-LH values
  //! 2. Both converge to stable values
  //! 3. Initial log-LH (before optimization) should be identical
  //! 4. Final log-LH difference should be bounded

  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::sync::{Arc, LazyLock};
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};

  pub static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  pub const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  pub fn gap_free_alignment() -> Result<Vec<FastaRecord>, Report> {
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

  pub fn setup_dense_only(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
  ) -> Result<Vec<Arc<RwLock<PartitionMarginalDense>>>, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
    let partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    initialize_marginal(graph, &partitions, aln)?;

    Ok(partitions)
  }

  pub fn setup_sparse_only(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
  ) -> Result<Vec<Arc<RwLock<PartitionMarginalSparse>>>, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
    let partitions = vec![Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(graph, &partitions, aln)?;
    update_marginal(graph, &partitions)?;

    Ok(partitions)
  }

  pub fn get_branch_lengths(graph: &GraphAncestral) -> Vec<f64> {
    graph
      .get_edges()
      .iter()
      .map(|edge| edge.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .collect()
  }
}
