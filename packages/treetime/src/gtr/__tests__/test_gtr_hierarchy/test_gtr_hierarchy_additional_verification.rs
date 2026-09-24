#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::AlphabetName;
  use crate::gtr::__tests__::rate_matrix::rate_matrix;
  use crate::gtr::get_gtr::{HKY85Params, JC69Params, K80Params, hky85, jc69, k80};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::array;

  #[test]
  fn test_gtr_jc69_q_symmetric() -> Result<(), Report> {
    let jc = jc69(JC69Params {
      mu: 1.0,
      alphabet: AlphabetName::Nuc,
    })?;

    let q = rate_matrix(&jc);
    assert_abs_diff_eq!(q, q.t(), epsilon = 1e-14);

    Ok(())
  }

  #[test]
  fn test_gtr_k80_kappa_1_q_symmetric() -> Result<(), Report> {
    let k = k80(K80Params {
      mu: 1.0,
      kappa: 1.0,
      alphabet: AlphabetName::Nuc,
    })?;

    let q = rate_matrix(&k);
    assert_abs_diff_eq!(q, q.t(), epsilon = 1e-14);

    Ok(())
  }

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

    let q1 = rate_matrix(&h1);
    let q2 = rate_matrix(&h2);
    assert_abs_diff_eq!(q1, q2, epsilon = 1e-14);

    assert!((h1.mu - h2.mu).abs() > 0.1);

    Ok(())
  }
}
