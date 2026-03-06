#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::AlphabetName;
  use crate::gtr::get_gtr::{HKY85Params, hky85};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::{Array2, array};

  // Branch Length Edge Cases

  /// At t = 0, no evolution occurs. The matrix exponential exp(Q*0) = I exactly.
  /// This is the base case for Felsenstein's pruning algorithm at leaf nodes.
  /// No tolerance needed - should be exact within floating-point representation.
  #[test]
  fn test_gtr_expqt_zero_branch() -> Result<(), Report> {
    let gtr = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: false,
    })?;

    let p = gtr.expQt(0.0);
    let identity = Array2::eye(4);

    // Should be exactly identity within floating-point precision
    assert_abs_diff_eq!(p, identity, epsilon = 1e-14);

    Ok(())
  }

  /// Very small branch lengths (t = 1e-10) test numerical stability near zero.
  /// The Taylor expansion exp(Qt) = I + Qt + (Qt)^2/2 + ... should give P ~ I.
  /// Risk: underflow in exp(lambda*t) for negative eigenvalues, or
  /// catastrophic cancellation in I + Qt when Qt is tiny.
  /// Verify: no NaN/Inf, P ~ I within tolerance.
  #[test]
  fn test_gtr_expqt_tiny_branch() -> Result<(), Report> {
    let gtr = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: false,
    })?;

    let t = 1e-10;
    let p = gtr.expQt(t);

    // Check no NaN or Inf
    assert!(
      !p.iter().any(|x| x.is_nan() || x.is_infinite()),
      "P contains NaN or Inf"
    );

    // P should be very close to identity
    let identity = Array2::eye(4);
    assert_abs_diff_eq!(p, identity, epsilon = 1e-8);

    Ok(())
  }

  /// For small t, the first-order Taylor approximation holds:
  ///
  ///   exp(Qt) ~ I + Qt + O(t^2)
  ///
  /// At t = 1e-6, higher-order terms contribute ~1e-12, negligible.
  /// This tests that the eigendecomposition method matches the direct approximation.
  #[test]
  fn test_gtr_expqt_small_branch_taylor() -> Result<(), Report> {
    let gtr = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(array![0.2, 0.3, 0.2, 0.3]),
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: false,
    })?;

    let t = 1e-6;
    let p = gtr.expQt(t);
    let q = gtr.Q();

    // First-order Taylor: P ~ I + mu*Q*t
    // Note: gtr.mu is already incorporated into eigenvalues, so we use mu*Q*t
    let taylor_approx = Array2::eye(4) + gtr.mu * t * &q;

    // At t=1e-6, error is O(t^2) ~ 1e-12
    assert_abs_diff_eq!(p, taylor_approx, epsilon = 1e-10);

    Ok(())
  }

  /// Large branch lengths (t = 100) should approach equilibrium.
  /// All eigenvalues except 0 have negative real parts, so:
  ///
  ///   exp(lambda_i * t) -> 0 as t -> inf (for lambda_i < 0)
  ///
  /// Only the stationary component survives: each column of P(t) -> pi.
  /// (In transposed convention, rows converge to pi.)
  ///
  /// Tolerance derivation (1e-6):
  /// - For t=100, mu=1, non-stationary terms decay as exp(-|lambda|*t)
  /// - GTR eigenvalues have |lambda| >= O(1), so exp(-100) ~ 1e-44
  /// - Remaining error from eigenvector reconstruction: V * diag(exp) * V_inv
  /// - Conservative bound 1e-6 accounts for ill-conditioned eigenvector matrices
  #[test]
  fn test_gtr_expqt_large_branch() -> Result<(), Report> {
    let gtr = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: false,
    })?;

    let t = 100.0;
    let p = gtr.expQt(t);

    // In transposed convention: P[i,j] -> pi[i] for all j
    // Each row of the equilibrium matrix is gtr.pi repeated across columns
    let expected = Array2::from_shape_fn((4, 4), |(i, _j)| gtr.pi[i]);

    // Tolerance 1e-6: see derivation in function docstring
    assert_abs_diff_eq!(p, expected, epsilon = 1e-6);

    Ok(())
  }

  /// At t = 1000, the system should have fully equilibrated.
  /// All rows of P(t) should equal pi within numerical precision.
  /// Risk: overflow in intermediate calculations if not using log-space
  /// or eigendecomposition properly.
  #[test]
  fn test_gtr_expqt_very_large_branch() -> Result<(), Report> {
    let gtr = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: false,
    })?;

    let t = 1000.0;
    let p = gtr.expQt(t);

    // Check no overflow (NaN or Inf)
    assert!(
      !p.iter().any(|x| x.is_nan() || x.is_infinite()),
      "P contains NaN or Inf at t={t}"
    );

    // Should be at equilibrium with high precision
    let expected = Array2::from_shape_fn((4, 4), |(i, _j)| gtr.pi[i]);
    assert_abs_diff_eq!(p, expected, epsilon = 1e-10);

    Ok(())
  }
}
