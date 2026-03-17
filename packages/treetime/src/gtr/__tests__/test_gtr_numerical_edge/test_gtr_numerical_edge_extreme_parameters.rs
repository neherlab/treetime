#[cfg(test)]
mod tests {
  use super::super::test_gtr_numerical_edge_support::assert_stochastic_matrix;
  use crate::alphabet::alphabet::AlphabetName;
  use crate::gtr::get_gtr::{HKY85Params, K80Params, hky85, k80};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::{Array2, array};

  // Extreme Parameter Tests

  /// kappa -> 0 means transitions (A<->G, C<->T) are strongly disfavored
  /// relative to transversions. The model remains valid but transitions
  /// become rare events. Tests that eigendecomposition handles near-singular
  /// rate matrices where some rates are much smaller than others.
  #[test]
  fn test_gtr_k80_kappa_near_zero() -> Result<(), Report> {
    let gtr = k80(K80Params {
      mu: 1.0,
      kappa: 0.01, // transitions 100x less likely than transversions
      alphabet: AlphabetName::Nuc,
    })?;

    let p = gtr.expQt(1.0);
    assert_stochastic_matrix(&p, "K80 kappa=0.01");

    Ok(())
  }

  /// kappa = 100 means transitions dominate transversions by 100:1.
  /// Biologically realistic for some organisms (e.g., vertebrate mitochondria).
  /// Tests numerical stability when rate matrix has widely varying magnitudes.
  /// The spectral gap may become very small, slowing convergence to equilibrium.
  #[test]
  fn test_gtr_k80_kappa_large() -> Result<(), Report> {
    let gtr = k80(K80Params {
      mu: 1.0,
      kappa: 100.0,
      alphabet: AlphabetName::Nuc,
    })?;

    let p = gtr.expQt(1.0);
    assert_stochastic_matrix(&p, "K80 kappa=100");

    // Transitions should be much more likely than transversions
    // A (index 0) <-> G (index 2) is a transition
    // A (index 0) <-> C (index 1) is a transversion
    // In P(t), P[dest, src] = prob of going from src to dest
    assert!(
      p[[2, 0]] > p[[1, 0]],
      "Expected P[G,A] > P[C,A] for large kappa: {} vs {}",
      p[[2, 0]],
      p[[1, 0]]
    );

    Ok(())
  }

  /// Highly skewed frequencies (pi = [0.97, 0.01, 0.01, 0.01]) occur in
  /// AT-rich or GC-rich genomes. Tests that:
  /// - Eigendecomposition handles near-singular matrices (small pi values)
  /// - sqrt(pi) in symmetrization doesn't cause numerical issues
  /// - Detailed balance still holds despite magnitude differences
  #[test]
  fn test_gtr_hky85_skewed_pi() -> Result<(), Report> {
    let pi = array![0.97, 0.01, 0.01, 0.01];

    let gtr = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(pi),
      alphabet: AlphabetName::Nuc,
    })?;

    let p = gtr.expQt(1.0);
    assert_stochastic_matrix(&p, "HKY85 skewed pi");

    // All entries bounded above
    assert!(
      !p.iter().any(|&x| x > 1.0 + 1e-14),
      "HKY85 skewed pi: matrix entry exceeds 1"
    );

    // Check detailed balance (key mathematical property)
    let q = gtr.Q();
    for i in 0..4 {
      for j in 0..4 {
        if i != j {
          let flux_ji = gtr.pi[j] * q[[i, j]];
          let flux_ij = gtr.pi[i] * q[[j, i]];
          assert_abs_diff_eq!(flux_ji, flux_ij, epsilon = 1e-10);
        }
      }
    }

    Ok(())
  }

  /// Nearly uniform pi (each ~ 0.25 + small perturbation) tests behavior
  /// near the symmetric JC69 limit. Small perturbations shouldn't cause
  /// large changes in the Q matrix or numerical instability.
  #[test]
  fn test_gtr_hky85_nearly_uniform_pi() -> Result<(), Report> {
    let pi = array![0.24, 0.26, 0.25, 0.25];

    let gtr = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(pi),
      alphabet: AlphabetName::Nuc,
    })?;

    let p = gtr.expQt(1.0);
    assert_stochastic_matrix(&p, "HKY85 nearly uniform pi");

    // Should be close to K80 with same kappa
    let k = k80(K80Params {
      mu: 1.0,
      kappa: 2.0,
      alphabet: AlphabetName::Nuc,
    })?;
    let p_k = k.expQt(1.0);

    // Verify continuity: small d_pi -> small d_P (Lipschitz-like bound).
    // Derivation: Q_ij ~ W_ij * pi_j, so d_Q_ij = W_ij * d_pi_j.
    // max(W_ij) = kappa = 2, max(d_pi) = 0.01, t = 1.
    // For P = exp(Qt), first-order d_P ~ t * d_Q, so max(d_P) ~ kappa * max(d_pi) = 0.02.
    // Observed ~0.007, well within bound.
    let diff = &p - &p_k;
    let max_diff = diff.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    let max_pi_diff = gtr.pi.iter().map(|x| (x - 0.25).abs()).fold(0.0_f64, f64::max);
    let kappa = 2.0;
    let bound = kappa * max_pi_diff; // 2.0 * 0.01 = 0.02

    assert!(
      max_diff < bound,
      "d_P exceeds Lipschitz bound: max_diff={max_diff:.6}, bound={bound:.4}, max_d_pi={max_pi_diff:.4}"
    );

    Ok(())
  }

  /// HKY85 must reduce to K80 when pi is exactly uniform.
  #[test]
  fn test_gtr_hky85_uniform_pi_matches_k80() -> Result<(), Report> {
    let h = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(array![0.25, 0.25, 0.25, 0.25]),
      alphabet: AlphabetName::Nuc,
    })?;
    let k = k80(K80Params {
      mu: 1.0,
      kappa: 2.0,
      alphabet: AlphabetName::Nuc,
    })?;

    assert_abs_diff_eq!(h.Q(), k.Q(), epsilon = 1e-14);
    assert_abs_diff_eq!(h.expQt(1.0), k.expQt(1.0), epsilon = 1e-14);

    Ok(())
  }

  /// mu = 0.001 gives very slow evolution (branch length in time >> in substitutions).
  /// exp(mu * Q * t) ~ I + mu*Q*t for moderate t.
  /// Tests that the implementation handles the slow-evolution regime without
  /// underflow in exp(small negative number).
  #[test]
  fn test_gtr_mu_very_small() -> Result<(), Report> {
    let gtr = hky85(HKY85Params {
      mu: 0.001,
      kappa: 2.0,
      pi: Some(array![0.25, 0.25, 0.25, 0.25]),
      alphabet: AlphabetName::Nuc,
    })?;

    let t = 1.0;
    let p = gtr.expQt(t);

    // With small mu, P should be close to identity
    for i in 0..4 {
      for j in 0..4 {
        let expected = if i == j { 1.0 } else { 0.0 };

        // mu*t = 0.001, so off-diagonal ~ 0.001, diagonal ~ 0.999
        assert!(
          (p[[i, j]] - expected).abs() < 0.01,
          "P[{i},{j}] = {} far from {}",
          p[[i, j]],
          expected
        );
      }
    }

    Ok(())
  }

  /// mu = 10.0 gives fast evolution (many substitutions per unit time).
  /// Even short branches approach equilibrium quickly.
  /// Tests that eigenvalue scaling doesn't cause overflow in exp(mu * lambda * t).
  #[test]
  fn test_gtr_mu_large() -> Result<(), Report> {
    let gtr = hky85(HKY85Params {
      mu: 10.0,
      kappa: 2.0,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      alphabet: AlphabetName::Nuc,
    })?;

    let t = 1.0;
    let p = gtr.expQt(t);

    // Check no overflow
    assert!(
      !p.iter().any(|x| x.is_nan() || x.is_infinite()),
      "P contains NaN or Inf"
    );

    // With large mu, should be near equilibrium even at t=1
    let expected = Array2::from_shape_fn((4, 4), |(i, _j)| gtr.pi[i]);
    assert_abs_diff_eq!(p, expected, epsilon = 1e-3);

    Ok(())
  }
}
