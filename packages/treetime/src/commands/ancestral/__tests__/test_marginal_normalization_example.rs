#[cfg(test)]
mod tests {
  use crate::commands::ancestral::__tests__::prop_generators::input::MarginalTestInput;
  use crate::commands::ancestral::__tests__::prop_marginal_support::tests::{run_dense_marginal, run_sparse_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use treetime_io::fasta::FastaRecord;
  use treetime_io::fasta::read_many_fasta_str;

  /// Build a fixed 4-taxon test input for marginal normalization verification.
  ///
  /// Tree topology: `((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01`
  /// Alignment: 8-character nucleotide sequences with shared prefix ACGTACG and all four
  /// nucleotides (T, A, G, C) at position 8, producing one variable site with maximum
  /// state diversity and seven fixed (invariant) sites.
  /// Model: JC69 (Jukes-Cantor 1969) with equal equilibrium frequencies (pi = 1/4).
  fn example_input() -> Result<MarginalTestInput, Report> {
    let alignment = read_many_fasta_str(
      "
>A
ACGTACGT
>B
ACGTACGA
>C
ACGTACGG
>D
ACGTACGC
",
      &*crate::test_utils::NUC_ALPHABET,
    )?;
    let gtr = jc69(JC69Params::default())?;
    Ok(MarginalTestInput {
      newick: "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;".to_owned(),
      alignment,
      gtr,
      n_taxa: 4,
      seq_len: 8,
    })
  }

  /// Validate that the example alignment has the expected shape: 4 sequences named
  /// A, B, C, D, each of length 8. Guards against accidental corruption of test fixtures.
  fn assert_example_alignment_shape(alignment: &[FastaRecord]) -> Result<(), Report> {
    let expected_names = ["A", "B", "C", "D"];
    for (index, record) in alignment.iter().enumerate() {
      let expected_name = expected_names[index];
      if record.seq_name != expected_name {
        return Err(eyre::eyre!(
          "Unexpected sequence name at index {index}: got {}, expected {expected_name}",
          record.seq_name
        ));
      }
      let expected_len = 8;
      let actual_len = record.seq.len();
      if actual_len != expected_len {
        return Err(eyre::eyre!(
          "Unexpected sequence length for {}: got {actual_len}, expected {expected_len}",
          record.seq_name
        ));
      }
    }
    Ok(())
  }

  /// Verify that Felsenstein's pruning algorithm (two-pass sum-product message passing on
  /// the tree) produces normalized marginal posterior distributions in the dense
  /// representation.
  ///
  /// The algorithm computes marginal posteriors P(x_k | D) at each node k via two passes:
  /// 1. Backward (leaves to root): partial likelihoods propagate upward (ingroup messages).
  /// 2. Forward (root to leaves): outgroup messages propagate downward.
  ///    The marginal posterior at each node is the product of ingroup and outgroup messages,
  ///    normalized to sum to 1 over all states at each alignment position.
  ///
  /// Checked invariants:
  /// - Log-likelihood is finite and non-positive (L = sum_s ln P(D_s) <= 0, since each
  ///   site likelihood P(D_s) <= 1).
  /// - Every row of every node profile matrix sums to 1.0 within 1e-8 (valid marginal
  ///   posterior distribution). Nodes without profiles (not yet populated) are skipped.
  /// - Every row of every edge outgroup message matrix sums to 1.0 within 1e-8
  ///   (normalized message from the rest of the tree toward the child subtree).
  /// - All values are finite and non-negative (within -1e-14 tolerance for numerical
  ///   noise from normalization).
  ///
  /// Example-based companion to `test_prop_marginal_normalization_dense`.
  /// Uses a fixed 4-taxon tree with JC69 model.
  #[test]
  fn test_marginal_normalization_example_dense() -> Result<(), Report> {
    let input = example_input()?;
    assert_example_alignment_shape(&input.alignment)?;
    let (log_lh, partitions) = run_dense_marginal(&input)?;
    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      if node_data.profile.dis.is_empty() {
        continue;
      }
      for row in node_data.profile.dis.rows() {
        let sum: f64 = row.sum();
        assert_abs_diff_eq!(1.0, sum, epsilon = 1e-8);
        for &value in row {
          assert!(value.is_finite(), "Non-finite value in dense node profile: {value}");
          assert!(value >= -1e-14, "Negative value in dense node profile: {value}");
        }
      }
    }
    for edge_data in partition.edges.values() {
      if edge_data.msg_to_child.dis.is_empty() {
        continue;
      }
      for row in edge_data.msg_to_child.dis.rows() {
        let sum: f64 = row.sum();
        assert_abs_diff_eq!(1.0, sum, epsilon = 1e-8);
        for &value in row {
          assert!(value.is_finite(), "Non-finite value in dense edge message: {value}");
          assert!(value >= -1e-14, "Negative value in dense edge message: {value}");
        }
      }
    }
    Ok(())
  }

  /// Verify that Felsenstein's pruning algorithm (two-pass sum-product message passing on
  /// the tree) produces normalized marginal posterior distributions in the sparse
  /// representation.
  ///
  /// The sparse representation uses Fitch parsimony compression to identify variable vs.
  /// fixed positions before running the marginal algorithm. Variable positions store
  /// individual probability vectors. Fixed positions are grouped by their conserved
  /// nucleotide character (the character that is invariant across the subtree at that
  /// node), each group sharing a single probability vector. Both components must be valid
  /// probability distributions after the two-pass algorithm (backward + forward).
  ///
  /// Checked invariants:
  /// - Log-likelihood is finite and non-positive (L = sum_s ln P(D_s) <= 0).
  /// - Node and edge profile `log_lh` fields are finite.
  /// - Every variable-position distribution sums to 1.0 within 1e-8.
  /// - Every fixed-character distribution sums to 1.0 within 1e-8.
  /// - All values are finite and non-negative (within -1e-14 tolerance for numerical
  ///   noise from normalization).
  ///
  /// Example-based companion to `test_prop_marginal_normalization_sparse`.
  /// Uses the same fixed 4-taxon tree and JC69 model as the dense variant.
  #[test]
  fn test_marginal_normalization_example_sparse() -> Result<(), Report> {
    let input = example_input()?;
    let (log_lh, partitions) = run_sparse_marginal(&input)?;
    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = partitions[0].read_arc();
    for node_data in partition.nodes.values() {
      let profile = &node_data.profile;
      assert!(
        profile.log_lh.is_finite(),
        "Sparse node profile log-lh is not finite: {}",
        profile.log_lh
      );
      for var_pos in profile.variable.values() {
        let sum: f64 = var_pos.dis.sum();
        assert_abs_diff_eq!(1.0, sum, epsilon = 1e-8);
        for &value in &var_pos.dis {
          assert!(
            value.is_finite(),
            "Non-finite value in sparse variable profile: {value}"
          );
          assert!(value >= -1e-14, "Negative value in sparse variable profile: {value}");
        }
      }
      for fixed_dis in profile.fixed.values() {
        let sum: f64 = fixed_dis.sum();
        assert_abs_diff_eq!(1.0, sum, epsilon = 1e-8);
        for &value in fixed_dis {
          assert!(value.is_finite(), "Non-finite value in sparse fixed profile: {value}");
          assert!(value >= -1e-14, "Negative value in sparse fixed profile: {value}");
        }
      }
    }
    for edge_data in partition.edges.values() {
      let profile = &edge_data.msg_to_child;
      assert!(
        profile.log_lh.is_finite(),
        "Sparse edge message log-lh is not finite: {}",
        profile.log_lh
      );
      for var_pos in profile.variable.values() {
        let sum: f64 = var_pos.dis.sum();
        assert_abs_diff_eq!(1.0, sum, epsilon = 1e-8);
        for &value in &var_pos.dis {
          assert!(
            value.is_finite(),
            "Non-finite value in sparse edge variable profile: {value}"
          );
          assert!(
            value >= -1e-14,
            "Negative value in sparse edge variable profile: {value}"
          );
        }
      }
      for fixed_dis in profile.fixed.values() {
        let sum: f64 = fixed_dis.sum();
        assert_abs_diff_eq!(1.0, sum, epsilon = 1e-8);
        for &value in fixed_dis {
          assert!(
            value.is_finite(),
            "Non-finite value in sparse edge fixed profile: {value}"
          );
          assert!(value >= -1e-14, "Negative value in sparse edge fixed profile: {value}");
        }
      }
    }
    Ok(())
  }
}
