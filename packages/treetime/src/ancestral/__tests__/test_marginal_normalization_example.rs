#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_generators::input::MarginalTestInput;
  use crate::ancestral::__tests__::prop_marginal_support::tests::{run_dense_marginal, run_sparse_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::{pretty_assert_array_finite, pretty_assert_array_nonneg};

  fn example_input() -> Result<MarginalTestInput, Report> {
    let alignment: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();
    let gtr = jc69(JC69Params::default())?;
    Ok(MarginalTestInput {
      newick: "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;".to_owned(),
      alignment,
      gtr,
      n_taxa: 4,
      seq_len: 8,
    })
  }

  fn assert_example_alignment_shape(alignment: &[AlignmentRecord]) -> Result<(), Report> {
    let expected_names = ["A", "B", "C", "D"];
    for (index, record) in alignment.iter().enumerate() {
      let expected_name = expected_names[index];
      if record.name != expected_name {
        return Err(eyre::eyre!(
          "Unexpected sequence name at index {index}: got {}, expected {expected_name}",
          record.name
        ));
      }
      let expected_len = 8;
      let actual_len = record.seq.len();
      if actual_len != expected_len {
        return Err(eyre::eyre!(
          "Unexpected sequence length for {}: got {actual_len}, expected {expected_len}",
          record.name
        ));
      }
    }
    Ok(())
  }

  #[test]
  fn test_marginal_normalization_example_dense() -> Result<(), Report> {
    let input = example_input()?;
    assert_example_alignment_shape(&input.alignment)?;
    let (log_lh, partitions) = run_dense_marginal(&input)?;
    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = &partitions;
    for node_data in partition.node_states.values() {
      if node_data.profile.dis.is_empty() {
        continue;
      }
      let dis = &node_data.profile.dis;
      pretty_assert_array_finite!(dis);
      pretty_assert_array_nonneg!(dis, epsilon = 1e-14);
      for row in dis.rows() {
        assert_abs_diff_eq!(1.0, row.sum(), epsilon = 1e-8);
      }
    }
    for edge_data in partition.edges.forward.values() {
      if edge_data.msg_to_child.dis.is_empty() {
        continue;
      }
      let dis = &edge_data.msg_to_child.dis;
      pretty_assert_array_finite!(dis);
      pretty_assert_array_nonneg!(dis, epsilon = 1e-14);
      for row in dis.rows() {
        assert_abs_diff_eq!(1.0, row.sum(), epsilon = 1e-8);
      }
    }
    Ok(())
  }

  #[test]
  fn test_marginal_normalization_example_sparse() -> Result<(), Report> {
    let input = example_input()?;
    let (log_lh, partitions) = run_sparse_marginal(&input)?;
    assert!(log_lh.is_finite(), "Log-likelihood is not finite: {log_lh}");
    assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");

    let partition = &partitions;
    for node_data in partition.node_states.values() {
      let profile = &node_data.profile;
      assert!(
        profile.log_lh.value().is_finite(),
        "Sparse node profile log-lh is not finite: {}",
        profile.log_lh.value()
      );
      for var_pos in profile.variable.values() {
        pretty_assert_array_finite!(var_pos.dis);
        pretty_assert_array_nonneg!(var_pos.dis, epsilon = 1e-14);
        assert_abs_diff_eq!(1.0, var_pos.dis.sum(), epsilon = 1e-8);
      }
      for fixed_dis in profile.fixed.values() {
        pretty_assert_array_finite!(fixed_dis);
        pretty_assert_array_nonneg!(fixed_dis, epsilon = 1e-14);
        assert_abs_diff_eq!(1.0, fixed_dis.sum(), epsilon = 1e-8);
      }
    }
    for edge_data in partition.edges.forward.values() {
      let profile = &edge_data.msg_to_child;
      assert!(
        profile.log_lh.value().is_finite(),
        "Sparse edge message log-lh is not finite: {}",
        profile.log_lh.value()
      );
      for var_pos in profile.variable.values() {
        pretty_assert_array_finite!(var_pos.dis);
        pretty_assert_array_nonneg!(var_pos.dis, epsilon = 1e-14);
        assert_abs_diff_eq!(1.0, var_pos.dis.sum(), epsilon = 1e-8);
      }
      for fixed_dis in profile.fixed.values() {
        pretty_assert_array_finite!(fixed_dis);
        pretty_assert_array_nonneg!(fixed_dis, epsilon = 1e-14);
        assert_abs_diff_eq!(1.0, fixed_dis.sum(), epsilon = 1e-8);
      }
    }
    Ok(())
  }
}
