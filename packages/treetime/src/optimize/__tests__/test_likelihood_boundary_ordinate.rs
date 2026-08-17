#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::likelihood::{branch_length_boundary_ordinate, evaluate_mixed_log_lh_only};
  use crate::partition::optimization_contribution::OptimizationContribution;
  use crate::partition::optimize_dense;
  use eyre::Report;
  use ndarray::array;
  use pretty_assertions::assert_eq;

  /// Single-site dense JC69 contribution whose per-site likelihood at `t = 0` is the coefficient row
  /// sum (every eigen-term carries `exp(eigval * 0) = 1`). A positive row sum gives a finite boundary
  /// log-likelihood; an all-zero row gives `L(0) = 0`, i.e. `ln L(0) = -inf`.
  fn dense_jc69_contribution(coefficient_row: [f64; 4]) -> OptimizationContribution {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![coefficient_row];
    OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr))
  }

  /// Finite boundary density: with no indels and a positive `L(0)`, the ordinate is the boundary
  /// neg-log peak-normalized against `max_log_lh`, i.e. `Some(max_log_lh - ln L(0))`. Oracle: the
  /// function contract's peak-normalization identity, with `ln L(0)` from the independent
  /// `evaluate_mixed_log_lh_only`.
  #[test]
  fn test_likelihood_boundary_ordinate_finite_returns_peak_normalized() -> Result<(), Report> {
    let contributions = [dense_jc69_contribution([0.0, 2.0, 0.0, 0.0])];
    let boundary_log_lh = evaluate_mixed_log_lh_only(&contributions, 0.0)?.value();
    let max_log_lh = 3.5;

    let actual = branch_length_boundary_ordinate(&contributions, 0, 0.0, max_log_lh)?;
    assert_eq!(Some(max_log_lh - boundary_log_lh), actual);
    Ok(())
  }

  /// Divergent boundary density from a degenerate site: `L(0) = 0` makes `ln L(0) = -inf`, so the
  /// ordinate is `None` (the caller must fit the power-law regime, not a straight line).
  #[test]
  fn test_likelihood_boundary_ordinate_divergent_site_returns_none() -> Result<(), Report> {
    let contributions = [dense_jc69_contribution([0.0, 0.0, 0.0, 0.0])];

    let actual = branch_length_boundary_ordinate(&contributions, 0, 0.0, 1.0)?;
    assert_eq!(None, actual);
    Ok(())
  }

  /// Indels short-circuit to `None` without evaluating the likelihood at `t = 0`: the same finite
  /// contribution that yields `Some` at `indel_count = 0` yields `None` once `indel_count > 0`,
  /// proving the short-circuit fires before the finiteness test.
  #[test]
  fn test_likelihood_boundary_ordinate_indel_short_circuits_to_none() -> Result<(), Report> {
    let contributions = [dense_jc69_contribution([0.0, 2.0, 0.0, 0.0])];

    let finite = branch_length_boundary_ordinate(&contributions, 0, 0.0, 3.5)?;
    assert!(finite.is_some(), "precondition: this contribution is finite at the boundary");

    let with_indel = branch_length_boundary_ordinate(&contributions, 1, 44.4, 3.5)?;
    assert_eq!(None, with_indel);
    Ok(())
  }
}
