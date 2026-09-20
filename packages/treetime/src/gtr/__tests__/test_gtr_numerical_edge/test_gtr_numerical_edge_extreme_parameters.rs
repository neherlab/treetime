#[cfg(test)]
mod tests {
  use super::super::test_gtr_numerical_edge_support::assert_stochastic_matrix;
  use crate::alphabet::alphabet::AlphabetName;
  use crate::gtr::get_gtr::{HKY85Params, K80Params, hky85, k80};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::{Array2, Axis, array};
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[test]
  fn test_gtr_k80_kappa_near_zero() -> Result<(), Report> {
    let gtr = k80(K80Params {
      mu: 1.0,
      kappa: 0.01,
      alphabet: AlphabetName::Nuc,
    })?;

    let p = gtr.expQt(1.0);
    assert_stochastic_matrix(&p, "K80 kappa=0.01");

    Ok(())
  }

  #[test]
  fn test_gtr_k80_kappa_large() -> Result<(), Report> {
    let gtr = k80(K80Params {
      mu: 1.0,
      kappa: 100.0,
      alphabet: AlphabetName::Nuc,
    })?;

    let p = gtr.expQt(1.0);
    assert_stochastic_matrix(&p, "K80 kappa=100");

    assert!(
      p[[2, 0]] > p[[1, 0]],
      "Expected P[G,A] > P[C,A] for large kappa: {} vs {}",
      p[[2, 0]],
      p[[1, 0]]
    );

    Ok(())
  }

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

    assert!(
      !p.iter().any(|&x| x > 1.0 + 1e-14),
      "HKY85 skewed pi: matrix entry exceeds 1"
    );

    let q = gtr.Q();
    let flux = &q * &gtr.pi.view().insert_axis(Axis(0));
    pretty_assert_abs_diff_eq!(flux, flux.t().to_owned(), epsilon = 1e-10);

    Ok(())
  }

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

    let k = k80(K80Params {
      mu: 1.0,
      kappa: 2.0,
      alphabet: AlphabetName::Nuc,
    })?;
    let p_k = k.expQt(1.0);

    let diff = &p - &p_k;
    let max_diff = diff.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    let max_pi_diff = gtr.pi.iter().map(|x| (x - 0.25).abs()).fold(0.0_f64, f64::max);
    let kappa = 2.0;
    let bound = kappa * max_pi_diff;

    assert!(
      max_diff < bound,
      "d_P exceeds Lipschitz bound: max_diff={max_diff:.6}, bound={bound:.4}, max_d_pi={max_pi_diff:.4}"
    );

    Ok(())
  }

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

    pretty_assert_abs_diff_eq!(p, Array2::eye(4), epsilon = 1e-2);

    Ok(())
  }

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

    assert!(
      !p.iter().any(|x| x.is_nan() || x.is_infinite()),
      "P contains NaN or Inf"
    );

    let expected = Array2::from_shape_fn((4, 4), |(i, _j)| gtr.pi[i]);
    assert_abs_diff_eq!(p, expected, epsilon = 1e-3);

    Ok(())
  }
}
