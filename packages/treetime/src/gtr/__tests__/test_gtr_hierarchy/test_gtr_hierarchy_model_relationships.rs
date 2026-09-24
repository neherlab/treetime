#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::__tests__::rate_matrix::rate_matrix;
  use crate::gtr::get_gtr::{F81Params, HKY85Params, JC69Params, K80Params, TN93Params, f81, hky85, jc69, k80, tn93};
  use crate::gtr::gtr::GTR;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::array;

  #[test]
  fn test_gtr_jc69_equals_k80_kappa_1() -> Result<(), Report> {
    let jc = jc69(JC69Params {
      mu: 1.0,
      alphabet: AlphabetName::Nuc,
    })?;

    let k = k80(K80Params {
      mu: 1.0,
      kappa: 1.0,
      alphabet: AlphabetName::Nuc,
    })?;

    let q_jc = rate_matrix(&jc);
    let q_k = rate_matrix(&k);
    assert_abs_diff_eq!(q_jc, q_k, epsilon = 1e-14);

    assert_abs_diff_eq!(jc.pi, k.pi, epsilon = 1e-14);

    Ok(())
  }

  #[test]
  fn test_gtr_jc69_equals_f81_uniform_pi() -> Result<(), Report> {
    let jc = jc69(JC69Params {
      mu: 1.0,
      alphabet: AlphabetName::Nuc,
    })?;

    let f = f81(F81Params {
      mu: 1.0,
      pi: Some(array![0.25, 0.25, 0.25, 0.25]),
      alphabet: AlphabetName::Nuc,
    })?;

    let q_jc = rate_matrix(&jc);
    let q_f = rate_matrix(&f);
    assert_abs_diff_eq!(q_jc, q_f, epsilon = 1e-14);

    assert_abs_diff_eq!(jc.pi, f.pi, epsilon = 1e-14);

    Ok(())
  }

  #[test]
  fn test_gtr_k80_equals_hky85_uniform_pi() -> Result<(), Report> {
    let kappa = 2.5;

    let k = k80(K80Params {
      mu: 1.0,
      kappa,
      alphabet: AlphabetName::Nuc,
    })?;

    let h = hky85(HKY85Params {
      mu: 1.0,
      kappa,
      pi: Some(array![0.25, 0.25, 0.25, 0.25]),
      alphabet: AlphabetName::Nuc,
    })?;

    let q_k = rate_matrix(&k);
    let q_h = rate_matrix(&h);
    assert_abs_diff_eq!(q_k, q_h, epsilon = 1e-14);

    assert_abs_diff_eq!(k.pi, h.pi, epsilon = 1e-14);

    Ok(())
  }

  #[test]
  fn test_gtr_f81_equals_hky85_kappa_1() -> Result<(), Report> {
    let pi = array![0.1, 0.2, 0.3, 0.4];

    let f = f81(F81Params {
      mu: 1.0,
      pi: Some(pi.clone()),
      alphabet: AlphabetName::Nuc,
    })?;

    let h = hky85(HKY85Params {
      mu: 1.0,
      kappa: 1.0,
      pi: Some(pi),
      alphabet: AlphabetName::Nuc,
    })?;

    let q_f = rate_matrix(&f);
    let q_h = rate_matrix(&h);
    assert_abs_diff_eq!(q_f, q_h, epsilon = 1e-14);

    assert_abs_diff_eq!(f.pi, h.pi, epsilon = 1e-14);

    Ok(())
  }

  #[test]
  fn test_gtr_hky85_equals_tn93_equal_transitions() -> Result<(), Report> {
    let pi = array![0.1, 0.2, 0.3, 0.4];

    let h = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(pi.clone()),
      alphabet: AlphabetName::Nuc,
    })?;

    let t = tn93(TN93Params {
      mu: 1.0,
      kappa1: 0.5,
      kappa2: 1.0,
      pi: Some(pi),
      alphabet: AlphabetName::Nuc,
    })?;

    let q_h = rate_matrix(&h);
    let q_t = rate_matrix(&t);
    assert_abs_diff_eq!(q_h, q_t, epsilon = 1e-14);

    assert_abs_diff_eq!(h.pi, t.pi, epsilon = 1e-14);

    Ok(())
  }

  #[test]
  fn test_gtr_tn93_equals_gtr_with_structured_w() -> Result<(), Report> {
    let pi = array![0.1, 0.2, 0.3, 0.4];
    let kappa1 = 0.5;
    let kappa2 = 1.5;

    let t = tn93(TN93Params {
      mu: 1.0,
      kappa1,
      kappa2,
      pi: Some(pi.clone()),
      alphabet: AlphabetName::Nuc,
    })?;

    #[rustfmt::skip]
    let w = array![
      [0.0,    kappa1, 1.0,    kappa1],
      [kappa1, 0.0,    kappa1, kappa2],
      [1.0,    kappa1, 0.0,    kappa1],
      [kappa1, kappa2, kappa1, 0.0   ]
    ];

    let alphabet = Alphabet::new(AlphabetName::Nuc).expect("Nuc alphabet should be valid");
    let n_states = alphabet.n_canonical();

    let g = GTR::builder().n_states(n_states).mu(1.0).W(w).pi(pi).build()?;

    let q_t = rate_matrix(&t);
    let q_g = rate_matrix(&g);
    assert_abs_diff_eq!(q_t, q_g, epsilon = 1e-14);

    assert_abs_diff_eq!(t.pi, g.pi, epsilon = 1e-14);

    Ok(())
  }
}
