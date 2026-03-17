#[cfg(test)]
pub mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::gtr::GTR;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::representation::payload::dense::DenseSeqDis;
  use crate::representation::payload::sparse::MarginalSparseSeqDistribution;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::sync::{Arc, LazyLock};
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  /// Assert that a dense marginal profile is numerically stable.
  pub fn assert_dense_profile_stable(profile: &DenseSeqDis, max_ulps: u32) {
    assert!(
      profile.log_lh.is_finite(),
      "Profile log_lh is not finite: {}",
      profile.log_lh
    );

    for (pos, row) in profile.dis.outer_iter().enumerate() {
      let sum: f64 = row.sum();
      assert_ulps_eq!(sum, 1.0, max_ulps = max_ulps);
      for (idx, &val) in row.iter().enumerate() {
        assert!(
          val.is_finite(),
          "Position {pos}, index {idx} has non-finite value: {val}"
        );
        assert!(val >= -1e-15, "Position {pos}, index {idx} has negative value: {val}");
      }
    }
  }

  /// Assert that a sparse marginal profile is numerically stable.
  pub fn assert_sparse_profile_stable(profile: &MarginalSparseSeqDistribution, max_ulps: u32) {
    assert!(
      profile.log_lh.is_finite(),
      "Profile log_lh is not finite: {}",
      profile.log_lh
    );

    for (pos, var_pos) in &profile.variable {
      let sum: f64 = var_pos.dis.sum();
      assert_ulps_eq!(sum, 1.0, max_ulps = max_ulps);
      for (idx, &val) in var_pos.dis.iter().enumerate() {
        assert!(
          val.is_finite(),
          "Variable position {pos}, index {idx} has non-finite value: {val}"
        );
        assert!(
          val >= -1e-15,
          "Variable position {pos}, index {idx} has negative value: {val}"
        );
      }
    }

    for (char_key, fixed_dis) in &profile.fixed {
      let sum: f64 = fixed_dis.sum();
      assert_ulps_eq!(sum, 1.0, max_ulps = max_ulps);
      for (idx, &val) in fixed_dis.iter().enumerate() {
        assert!(
          val.is_finite(),
          "Fixed distribution for char {char_key:?}, index {idx} has non-finite value: {val}"
        );
        assert!(
          val >= -1e-15,
          "Fixed distribution for char {char_key:?}, index {idx} has negative value: {val}"
        );
      }
    }
  }

  /// Run dense marginal reconstruction with a custom GTR model and return the log-likelihood and partition array.
  pub fn run_dense_marginal_with_partitions(
    newick: &str,
    aln_str: &str,
    gtr: GTR,
  ) -> Result<(f64, [Arc<RwLock<PartitionMarginalDense>>; 1]), Report> {
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    let log_lh = initialize_marginal(&graph, &partitions, &aln)?;
    Ok((log_lh, partitions))
  }

  /// Run sparse marginal reconstruction with a custom GTR model and return the log-likelihood and partition array.
  pub fn run_sparse_marginal_with_partitions(
    newick: &str,
    aln_str: &str,
    gtr: GTR,
  ) -> Result<(f64, [Arc<RwLock<PartitionMarginalSparse>>; 1]), Report> {
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions, &aln)?;
    let log_lh = update_marginal(&graph, &partitions)?;
    Ok((log_lh, partitions))
  }
}
