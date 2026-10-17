#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{
    F81Params, HKY85Params, JC69Params, Jtt92Params, K80Params, T92Params, TN93Params, f81, hky85, jc69, jtt92, k80,
    t92, tn93,
  };
  use crate::gtr::gtr::{GTR, Mu};
  use approx::{assert_abs_diff_eq, assert_ulps_eq};
  use eyre::Report;
  use helpers::normalized;
  use ndarray::{Array1, Array2, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::jc69(          jc69(JC69Params::default()).unwrap(),                                                  array![0.25, 0.25, 0.25, 0.25])]
  #[case::f81_uniform(   f81(F81Params::default()).unwrap(),                                                    array![0.25, 0.25, 0.25, 0.25])]
  #[case::f81_nonuniform(f81(F81Params { pi: Some(array![0.1, 0.2, 0.3, 0.4]), ..F81Params::default() }).unwrap(), array![0.1, 0.2, 0.3, 0.4])]
  #[case::k80_kappa_one( k80(K80Params { kappa: 1.0, ..K80Params::default() }).unwrap(),                        array![0.25, 0.25, 0.25, 0.25])]
  #[case::tn93_default(  tn93(TN93Params::default()).unwrap(),                                                  array![0.25, 0.25, 0.25, 0.25])]
  #[case::binary(        normalized(array![0.3, 0.7], None).unwrap(),                                           array![0.3, 0.7])]
  #[case::three_states(  normalized(array![0.2, 0.3, 0.5], Some(Array2::from_elem((3, 3), 2.0))).unwrap(),      array![0.2, 0.3, 0.5])]
  #[trace]
  fn test_gtr_construction_equal_exchangeabilities_certify_f81_decay_rate(
    #[case] gtr: GTR,
    #[case] pi: Array1<f64>,
  ) {
    let f81_decay_rate = 1.0 / (1.0 - pi.dot(&pi));
    assert_abs_diff_eq!(f81_decay_rate, gtr.certified_decay_rate().unwrap(), epsilon = 1e-14);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::k80_default(  k80(K80Params::default()).unwrap())]
  #[case::hky85_default(hky85(HKY85Params::default()).unwrap())]
  #[case::t92_default(  t92(T92Params::default()).unwrap())]
  #[case::tn93_kappa(   tn93(TN93Params { kappa1: 2.0, ..TN93Params::default() }).unwrap())]
  #[case::jtt92(        jtt92(Jtt92Params::default()).unwrap())]
  #[case::one_unequal(  normalized(array![0.2, 0.3, 0.5], Some(array![[0.0, 1.0, 1.0], [1.0, 0.0, 1.5], [1.0, 1.5, 0.0]])).unwrap())]
  #[case::one_zero(     normalized(array![0.2, 0.3, 0.5], Some(array![[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [1.0, 1.0, 0.0]])).unwrap())]
  #[trace]
  fn test_gtr_construction_unequal_exchangeabilities_are_not_certified(#[case] gtr: GTR) {
    assert_eq!(None, gtr.certified_decay_rate());
  }

  #[test]
  fn test_gtr_construction_asymmetric_input_certified_after_symmetrization() -> Result<(), Report> {
    let gtr = normalized(
      array![0.25, 0.25, 0.25, 0.25],
      Some(array![
        [0.0, 0.5, 1.5, 1.0],
        [1.5, 0.0, 1.0, 0.5],
        [0.5, 1.0, 0.0, 1.5],
        [1.0, 1.5, 0.5, 0.0]
      ]),
    )?;
    let jc69_decay_rate = 4.0 / 3.0;
    assert_abs_diff_eq!(jc69_decay_rate, gtr.certified_decay_rate().unwrap(), epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_gtr_construction_decay_rate_matches_nonzero_eigenvalue() -> Result<(), Report> {
    let gtr = f81(F81Params {
      mu: 0.7,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      ..F81Params::default()
    })?;
    let decay_rate = gtr.certified_decay_rate().unwrap();
    let rate_eigvals = gtr.eigvals() * gtr.mu();
    let expected = array![-decay_rate, -decay_rate, -decay_rate, 0.0];
    assert_abs_diff_eq!(expected, rate_eigvals, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_gtr_construction_decay_rate_follows_mu_updates() -> Result<(), Report> {
    let mut gtr = jc69(JC69Params::default())?;
    gtr.set_mu(2.0)?;
    let jc69_decay_rate_at_mu_2 = 8.0 / 3.0;
    assert_abs_diff_eq!(
      jc69_decay_rate_at_mu_2,
      gtr.certified_decay_rate().unwrap(),
      epsilon = 1e-14
    );
    Ok(())
  }

  #[test]
  fn test_gtr_construction_decay_rate_from_matrix_scale() -> Result<(), Report> {
    let gtr = GTR::builder()
      .n_states(3)
      .mu(Mu::FromMatrix(1.5))
      .W(Array2::from_elem((3, 3), 2.0))
      .pi(array![0.2, 0.3, 0.5])
      .build()?;
    let mu_times_exchangeability = 1.5 * 2.0;
    assert_abs_diff_eq!(
      mu_times_exchangeability,
      gtr.certified_decay_rate().unwrap(),
      epsilon = 1e-14
    );
    Ok(())
  }

  #[test]
  fn test_gtr_construction_set_mu_keeps_eigensystem() -> Result<(), Report> {
    let original = hky85(HKY85Params::default())?;
    let mut updated = original.clone();
    updated.set_mu(3.0)?;
    assert_ulps_eq!(3.0, updated.mu(), max_ulps = 0);
    assert_eq!(original.eigvals(), updated.eigvals());
    assert_eq!(original.v(), updated.v());
    assert_eq!(original.v_inv(), updated.v_inv());
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::negative(-1.0,          "Substitution rate (mu) must be finite and non-negative, got -1")]
  #[case::nan(     f64::NAN,      "Substitution rate (mu) must be finite and non-negative, got NaN")]
  #[case::infinite(f64::INFINITY, "Substitution rate (mu) must be finite and non-negative, got inf")]
  #[trace]
  fn test_gtr_construction_rejects_invalid_mu(#[case] mu: f64, #[case] expected: &str) -> Result<(), Report> {
    assert_error!(
      GTR::builder().n_states(2).mu(Mu::Normalized(mu)).pi(array![0.5, 0.5]).build(),
      expected
    );
    let mut gtr = jc69(JC69Params::default())?;
    assert_error!(gtr.set_mu(mu), expected);
    Ok(())
  }

  #[test]
  fn test_gtr_construction_rejects_overflowing_resolved_mu() {
    assert_error!(
      GTR::builder()
        .n_states(2)
        .mu(Mu::FromMatrix(f64::MAX))
        .W(Array2::from_elem((2, 2), 4.0))
        .pi(array![0.5, 0.5])
        .build(),
      "Substitution rate (mu) must be finite and non-negative, got inf"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::negative(array![-0.5, 1.5],     "Equilibrium frequencies (pi) must be finite and positive with a finite sum, got [-0.5, 1.5]")]
  #[case::nan(     array![f64::NAN, 1.0], "Equilibrium frequencies (pi) must be finite and positive with a finite sum, got [NaN, 1]")]
  #[case::one_zero(array![0.0, 1.0],      "Equilibrium frequencies (pi) must be finite and positive with a finite sum, got [0, 1]")]
  #[case::all_zero(array![0.0, 0.0],      "Equilibrium frequencies (pi) must be finite and positive with a finite sum, got [0, 0]")]
  #[trace]
  fn test_gtr_construction_rejects_invalid_pi(#[case] pi: Array1<f64>, #[case] expected: &str) {
    assert_error!(normalized(pi, None), expected);
  }

  #[test]
  fn test_gtr_construction_rejects_pi_with_infinite_sum() {
    let pi = array![f64::MAX, f64::MAX];
    let expected = format!("Equilibrium frequencies (pi) must be finite and positive with a finite sum, got {pi}");
    assert_error!(normalized(pi, None), expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::negative(array![[0.0, -1.0], [-1.0, 0.0]],          "Substitution matrix (W) must contain finite non-negative entries, got [[0, -1],\n [-1, 0]]")]
  #[case::nan(     array![[0.0, f64::NAN], [f64::NAN, 0.0]], "Substitution matrix (W) must contain finite non-negative entries, got [[0, NaN],\n [NaN, 0]]")]
  #[case::zero(    array![[0.0, 0.0], [0.0, 0.0]],            "Average substitution rate of the substitution matrix (W) under the equilibrium frequencies (pi) must be positive, got 0")]
  #[trace]
  fn test_gtr_construction_rejects_invalid_w(#[case] W: Array2<f64>, #[case] expected: &str) {
    assert_error!(normalized(array![0.5, 0.5], Some(W)), expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::negative(array![1.0, -0.5],       "Site rates must be finite and non-negative, got -0.5 at site 1")]
  #[case::nan(     array![f64::NAN, 1.0],   "Site rates must be finite and non-negative, got NaN at site 0")]
  #[case::infinite(array![1.0, f64::INFINITY], "Site rates must be finite and non-negative, got inf at site 1")]
  #[trace]
  fn test_gtr_construction_rejects_invalid_site_rates(#[case] site_rates: Array1<f64>, #[case] expected: &str) {
    assert_error!(
      GTR::builder()
        .n_states(2)
        .mu(Mu::Normalized(1.0))
        .pi(array![0.5, 0.5])
        .site_rates(site_rates)
        .build(),
      expected
    );
  }

  #[test]
  fn test_gtr_construction_accepts_zero_site_rate() -> Result<(), Report> {
    let site_rates = array![0.0, 1.0];
    let gtr = GTR::builder()
      .n_states(2)
      .mu(Mu::Normalized(1.0))
      .pi(array![0.5, 0.5])
      .site_rates(site_rates.clone())
      .build()?;
    assert_eq!(Some(&site_rates), gtr.site_rates());
    Ok(())
  }

  mod helpers {
    use crate::gtr::gtr::{GTR, Mu};
    use eyre::Report;
    use ndarray::{Array1, Array2};

    pub(super) fn normalized(pi: Array1<f64>, W: Option<Array2<f64>>) -> Result<GTR, Report> {
      GTR::builder()
        .n_states(pi.len())
        .mu(Mu::Normalized(1.0))
        .maybe_W(W)
        .pi(pi)
        .build()
    }
  }
}
