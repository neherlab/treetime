#[cfg(test)]
mod tests {
  use crate::coalescent::skyline::{marginal_log_tc_variances, skyline_hessian};
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;

  #[test]
  fn test_skyline_marginal_variance_matches_2x2_inverse() -> Result<(), Report> {
    let (a, b, s) = (5.0_f64, 3.0_f64, 1.0_f64);
    let det = a * b - s.powi(2);

    let hessian = skyline_hessian(&[0.0, 0.0], &[a - s, b - s], s)?;
    let variances = marginal_log_tc_variances(&hessian)?;

    pretty_assert_ulps_eq!(b / det, variances[0], max_ulps = 8);
    pretty_assert_ulps_eq!(a / det, variances[1], max_ulps = 8);

    Ok(())
  }

  #[test]
  fn test_skyline_marginal_variance_exceeds_diagonal_only_under_coupling() -> Result<(), Report> {
    let (a, b, s) = (5.0, 3.0, 1.0);

    let hessian = skyline_hessian(&[0.0, 0.0], &[a - s, b - s], s)?;
    let variances = marginal_log_tc_variances(&hessian)?;

    assert!(
      variances[0] > 1.0 / a,
      "marginal variance {} must exceed the diagonal-only {}",
      variances[0],
      1.0 / a
    );
    assert!(
      variances[1] > 1.0 / b,
      "marginal variance {} must exceed the diagonal-only {}",
      variances[1],
      1.0 / b
    );

    Ok(())
  }
}
