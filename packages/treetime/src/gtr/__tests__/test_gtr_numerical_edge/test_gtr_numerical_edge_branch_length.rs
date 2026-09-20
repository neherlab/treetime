#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::AlphabetName;
  use crate::gtr::get_gtr::{HKY85Params, hky85};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::{Array2, array};

  #[test]
  fn test_gtr_expqt_zero_branch() -> Result<(), Report> {
    let gtr = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      alphabet: AlphabetName::Nuc,
    })?;

    let p = gtr.expQt(0.0);
    let identity = Array2::eye(4);

    assert_abs_diff_eq!(p, identity, epsilon = 1e-14);

    Ok(())
  }

  #[test]
  fn test_gtr_expqt_tiny_branch() -> Result<(), Report> {
    let gtr = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      alphabet: AlphabetName::Nuc,
    })?;

    let t = 1e-10;
    let p = gtr.expQt(t);

    assert!(
      !p.iter().any(|x| x.is_nan() || x.is_infinite()),
      "P contains NaN or Inf"
    );

    let identity = Array2::eye(4);
    assert_abs_diff_eq!(p, identity, epsilon = 1e-8);

    Ok(())
  }

  #[test]
  fn test_gtr_expqt_small_branch_taylor() -> Result<(), Report> {
    let gtr = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(array![0.2, 0.3, 0.2, 0.3]),
      alphabet: AlphabetName::Nuc,
    })?;

    let t = 1e-6;
    let p = gtr.expQt(t);
    let q = gtr.Q();

    let taylor_approx = Array2::eye(4) + gtr.mu * t * &q;

    assert_abs_diff_eq!(p, taylor_approx, epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_gtr_expqt_large_branch() -> Result<(), Report> {
    let gtr = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      alphabet: AlphabetName::Nuc,
    })?;

    let t = 100.0;
    let p = gtr.expQt(t);

    let expected = Array2::from_shape_fn((4, 4), |(i, _j)| gtr.pi[i]);

    assert_abs_diff_eq!(p, expected, epsilon = 1e-6);

    Ok(())
  }

  #[test]
  fn test_gtr_expqt_very_large_branch() -> Result<(), Report> {
    let gtr = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      alphabet: AlphabetName::Nuc,
    })?;

    let t = 1000.0;
    let p = gtr.expQt(t);

    assert!(
      !p.iter().any(|x| x.is_nan() || x.is_infinite()),
      "P contains NaN or Inf at t={t}"
    );

    let expected = Array2::from_shape_fn((4, 4), |(i, _j)| gtr.pi[i]);
    assert_abs_diff_eq!(p, expected, epsilon = 1e-10);

    Ok(())
  }
}
