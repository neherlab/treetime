#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::JC69Params;
  use crate::gtr::get_gtr::jc69;
  use crate::test_utils::{run_dense_marginal_with_newick, run_sparse_marginal_with_newick};
  use eyre::Report;
  use itertools::Itertools;

  /// Felsenstein site likelihood increases monotonically with branch length for mismatched leaf observations.
  #[test]
  fn test_likelihood_monotonic_increase_mismatched_sequences() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let branch_lengths = (1..=20).map(|i| i as f64 * 0.1).collect_vec();
    let aln = ">A\nA\n>B\nT\n";

    let log_lhs: Vec<f64> = branch_lengths
      .iter()
      .map(|t| {
        let newick = format!("(A:{t},B:{t})root;");
        run_dense_marginal_with_newick(&newick, aln, gtr.clone())
      })
      .collect::<Result<Vec<_>, _>>()?;

    let violation = branch_lengths
      .iter()
      .zip(log_lhs.iter())
      .tuple_windows()
      .find(|((_, prev_log_lh), (_, log_lh))| **log_lh < **prev_log_lh - 1e-10);

    if let Some(((t, prev_log_lh), (_, log_lh))) = violation {
      panic!(
        "Likelihood should increase monotonically for mismatched sequences. At t={t}: log_lh={log_lh}, prev={prev_log_lh}"
      );
    }

    let equilibrium_log_lh = (0.25 * 0.25_f64).ln();
    let final_log_lh = *log_lhs.last().unwrap_or(&f64::NEG_INFINITY);
    assert!(
      (final_log_lh - equilibrium_log_lh).abs() < 0.1,
      "At t=2.0, should be close to equilibrium. actual={final_log_lh}, equilibrium={equilibrium_log_lh}"
    );

    Ok(())
  }

  /// Felsenstein site likelihood decreases monotonically with branch length for matched leaf observations.
  #[test]
  fn test_likelihood_maximized_near_zero_for_matched_sequences() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let branch_lengths = (1..=20).map(|i| i as f64 * 0.1).collect_vec();

    let log_lhs: Vec<f64> = branch_lengths
      .iter()
      .map(|t| {
        let newick = format!("(A:{t},B:{t})root;");
        let aln = ">A\nA\n>B\nA\n";
        run_dense_marginal_with_newick(&newick, aln, gtr.clone())
      })
      .collect::<Result<Vec<_>, _>>()?;

    let violation = branch_lengths
      .iter()
      .zip(log_lhs.iter())
      .tuple_windows()
      .find(|((_, prev_log_lh), (_, log_lh))| **log_lh > **prev_log_lh + 1e-10);

    if let Some(((t, prev_log_lh), (_, log_lh))) = violation {
      panic!(
        "Likelihood should decrease as branch length increases for matched sequences. At t={t}: log_lh={log_lh}, prev={prev_log_lh}"
      );
    }

    Ok(())
  }

  /// Numerical stability: log-likelihood remains finite and non-positive across 4 orders of magnitude in branch length.
  #[test]
  fn test_likelihood_finite_across_branch_length_range() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let branch_lengths = [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0];

    let results: Vec<(f64, f64)> = branch_lengths
      .iter()
      .map(|t| {
        let newick = format!("(A:{t},B:{t})root;");
        let aln = ">A\nACGT\n>B\nTGCA\n";
        let log_lh = run_dense_marginal_with_newick(&newick, aln, gtr.clone())?;
        Ok::<_, Report>((*t, log_lh))
      })
      .collect::<Result<Vec<_>, _>>()?;

    let violation = results
      .into_iter()
      .find(|(_, log_lh)| !log_lh.is_finite() || *log_lh > 0.0);

    if let Some((t, log_lh)) = violation {
      panic!("Log-likelihood should stay finite and non-positive. At t={t}: log_lh={log_lh}");
    }

    Ok(())
  }

  /// Monotonic decrease of Felsenstein likelihood for identical sequences on a three-taxon tree with an internal node.
  #[test]
  fn test_likelihood_monotonicity_three_taxon_tree() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let branch_lengths = (1..=15).map(|i| i as f64 * 0.1).collect_vec();

    let log_lhs: Vec<f64> = branch_lengths
      .iter()
      .map(|t| {
        let newick = format!("((A:{t},B:{t})AB:{t},C:{t})root;");
        let aln = ">A\nAAAA\n>B\nAAAA\n>C\nAAAA\n";
        run_dense_marginal_with_newick(&newick, aln, gtr.clone())
      })
      .collect::<Result<Vec<_>, _>>()?;

    let violation = branch_lengths
      .iter()
      .zip(log_lhs.iter())
      .tuple_windows()
      .find(|((_, prev_log_lh), (_, log_lh))| **log_lh > **prev_log_lh + 1e-10);

    if let Some(((t, prev_log_lh), (_, log_lh))) = violation {
      panic!(
        "Likelihood should decrease for identical sequences as t increases. At t={t}: log_lh={log_lh}, prev={prev_log_lh}"
      );
    }

    Ok(())
  }

  /// Monotonic decrease of Felsenstein likelihood for identical sequences using sparse representation.
  #[test]
  fn test_likelihood_monotonicity_sparse_partition() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let branch_lengths = (1..=15).map(|i| i as f64 * 0.1).collect_vec();

    let log_lhs: Vec<f64> = branch_lengths
      .iter()
      .map(|t| {
        let newick = format!("(A:{t},B:{t})root;");
        let aln = ">A\nGGGG\n>B\nGGGG\n";
        run_sparse_marginal_with_newick(&newick, aln, gtr.clone())
      })
      .collect::<Result<Vec<_>, _>>()?;

    let violation = branch_lengths
      .iter()
      .zip(log_lhs.iter())
      .tuple_windows()
      .find(|((_, prev_log_lh), (_, log_lh))| **log_lh > **prev_log_lh + 1e-10);

    if let Some(((t, prev_log_lh), (_, log_lh))) = violation {
      panic!(
        "Sparse likelihood should decrease for identical sequences. At t={t}: log_lh={log_lh}, prev={prev_log_lh}"
      );
    }

    Ok(())
  }
}
