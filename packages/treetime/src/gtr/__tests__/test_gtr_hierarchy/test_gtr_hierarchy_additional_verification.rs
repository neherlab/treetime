#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::AlphabetName;
  use crate::gtr::get_gtr::{HKY85Params, JC69Params, K80Params, hky85, jc69, k80};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::array;

  // Additional Hierarchy Verification Tests

  /// Verify that JC69 produces a symmetric Q matrix (since pi is uniform and W is symmetric).
  #[test]
  fn test_gtr_jc69_q_symmetric() -> Result<(), Report> {
    let jc = jc69(JC69Params {
      mu: 1.0,
      alphabet: AlphabetName::Nuc,
    })?;

    let q = jc.Q();
    assert_abs_diff_eq!(q, q.t(), epsilon = 1e-14);

    Ok(())
  }

  /// K80 with kappa=1 should have symmetric Q (same as JC69).
  #[test]
  fn test_gtr_k80_kappa_1_q_symmetric() -> Result<(), Report> {
    let k = k80(K80Params {
      mu: 1.0,
      kappa: 1.0,
      alphabet: AlphabetName::Nuc,
    })?;

    let q = k.Q();
    assert_abs_diff_eq!(q, q.t(), epsilon = 1e-14);

    Ok(())
  }

  /// Verify that varying mu doesn't change the relative rates in Q.
  /// Q is normalized, so changing mu scales the eigenvalues but not Q itself.
  #[test]
  fn test_gtr_mu_does_not_affect_q_shape() -> Result<(), Report> {
    let h1 = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      alphabet: AlphabetName::Nuc,
    })?;

    let h2 = hky85(HKY85Params {
      mu: 5.0,
      kappa: 2.0,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      alphabet: AlphabetName::Nuc,
    })?;

    // Q matrices should be identical (mu only affects eigenvalues)
    let q1 = h1.Q();
    let q2 = h2.Q();
    assert_abs_diff_eq!(q1, q2, epsilon = 1e-14);

    // But mu should differ
    assert!((h1.mu - h2.mu).abs() > 0.1);

    Ok(())
  }
}
