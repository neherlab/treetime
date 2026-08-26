#[cfg(test)]
mod tests {
  use crate::coalescent::skyline::{marginal_log_tc_variances, skyline_hessian};
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;

  // The multi-segment marginal variance inverts the full penalized tridiagonal Hessian
  // and takes the diagonal, so a diagonal-only shortcut would pass any test that only
  // checks finiteness or bracketing. These two tests pin the exact inverse against a
  // hand-derived 2x2 oracle and the covariance contribution the off-diagonal coupling
  // adds, distinguishing the true implementation from the diagonal-only alternative.
  //
  // Construction: `skyline_hessian` builds `diag_i = I_i * e^{-z_i} + stiffness` and
  // off-diagonal `-stiffness`. With `z = [0, 0]` the data curvature is `I_i`, so choosing
  // `I = [a - s, b - s]` and `stiffness = s` yields exactly `H = [[a, -s], [-s, b]]`.

  #[test]
  fn test_skyline_marginal_variance_matches_2x2_inverse() -> Result<(), Report> {
    // Oracle: closed-form 2x2 inverse. For H = [[a, -s], [-s, b]], det = a*b - s^2 and
    // diag(H^{-1}) = [b/det, a/det].
    let (a, b, s) = (5.0, 3.0, 1.0);
    let det = a * b - s * s;

    let hessian = skyline_hessian(&[0.0, 0.0], &[a - s, b - s], s)?;
    let variances = marginal_log_tc_variances(&hessian)?;

    pretty_assert_ulps_eq!(b / det, variances[0], max_ulps = 8);
    pretty_assert_ulps_eq!(a / det, variances[1], max_ulps = 8);

    Ok(())
  }

  #[test]
  fn test_skyline_marginal_variance_exceeds_diagonal_only_under_coupling() -> Result<(), Report> {
    // Metamorphic: with s > 0 the segments are coupled, so each marginal variance is
    // strictly larger than the diagonal-only reciprocal `1/H_ii` it would carry in
    // isolation: b/(a*b - s^2) > 1/a  and  a/(a*b - s^2) > 1/b, since a*b > a*b - s^2.
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
