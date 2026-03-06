#[cfg(test)]
mod tests {
  use super::super::test_gtr_numerical_edge_support::assert_stochastic_matrix;
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::get_gtr::{HKY85Params, K80Params, hky85, k80};
  use crate::gtr::gtr::{GTR, GTRParams};
  use eyre::Report;
  use ndarray::array;
  use rstest::rstest;

  // Parameterized Edge Case Tests

  /// Test expQt at various branch lengths for stability.
  #[rstest]
  #[case::zero(0.0)]
  #[case::tiny(1e-15)]
  #[case::small(1e-10)]
  #[case::medium(1e-5)]
  #[case::one(1.0)]
  #[case::large(100.0)]
  #[case::very_large(1000.0)]
  fn test_gtr_expqt_branch_length_range(#[case] t: f64) -> Result<(), Report> {
    let gtr = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: false,
    })?;

    let p = gtr.expQt(t);

    // No NaN or Inf
    assert!(
      !p.iter().any(|x| x.is_nan() || x.is_infinite()),
      "P(t={t}) contains NaN or Inf"
    );
    assert_stochastic_matrix(&p, &format!("expQt t={t}"));

    Ok(())
  }

  /// Test various kappa values for K80 model stability.
  #[rstest]
  #[case::very_small(0.001)]
  #[case::small(0.1)]
  #[case::one(1.0)]
  #[case::moderate(5.0)]
  #[case::large(50.0)]
  #[case::very_large(500.0)]
  fn test_gtr_k80_kappa_range(#[case] kappa: f64) -> Result<(), Report> {
    let gtr = k80(K80Params {
      mu: 1.0,
      kappa,
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: false,
    })?;

    let p = gtr.expQt(1.0);
    assert!(!p.iter().any(|x| x.is_nan()), "P contains NaN for kappa={kappa}");
    assert_stochastic_matrix(&p, &format!("K80 kappa={kappa}"));

    Ok(())
  }

  /// Test custom GTR with extreme but valid W values.
  #[test]
  fn test_gtr_extreme_w_values() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;

    // W with very different magnitudes
    #[rustfmt::skip]
    let w = array![
      [0.0,   0.01, 100.0, 0.1  ],
      [0.01,  0.0,  0.5,   50.0 ],
      [100.0, 0.5,  0.0,   0.01 ],
      [0.1,   50.0, 0.01,  0.0  ]
    ];

    let pi = array![0.25, 0.25, 0.25, 0.25];
    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: Some(w),
      pi,
    })?;

    let p = gtr.expQt(1.0);
    assert!(
      !p.iter().any(|x| x.is_nan() || x.is_infinite()),
      "P contains NaN or Inf"
    );
    assert_stochastic_matrix(&p, "extreme W values");

    Ok(())
  }
}
