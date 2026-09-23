#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {

  use crate::gtr::get_gtr::{JC69Params, jc69};

  use crate::optimize::indel::poisson_indel_log_lh;
  use crate::optimize::likelihood::{evaluate_mixed, evaluate_mixed_log_lh_only};

  use crate::optimize::method_newton::{chain_rule_log, chain_rule_sqrt};

  use crate::partition::optimize;
  use crate::partition::optimize::contribution::OptimizationContribution;

  use approx::assert_abs_diff_eq;

  use ndarray::array;
  use rstest::rstest;

  #[test]
  fn test_optimize_method_chain_rule_at_zero() {
    let (ds, d2s) = chain_rule_sqrt(0.0, 100.0, -500.0);
    assert_abs_diff_eq!(ds, 0.0, epsilon = 1e-15);
    assert_abs_diff_eq!(d2s, 200.0, epsilon = 1e-15);
  }

  #[test]
  fn test_optimize_method_chain_rule_analytical() {
    let s = 0.3;
    let (ds, d2s) = chain_rule_sqrt(s, 10.0, -100.0);
    assert_abs_diff_eq!(ds, 6.0, epsilon = 1e-14);
    assert_abs_diff_eq!(d2s, -16.0, epsilon = 1e-14);
  }

  #[test]
  fn test_optimize_method_chain_rule_log_analytical() {
    let t = 0.09;
    let (du, d2u) = chain_rule_log(t, 10.0, -100.0);
    assert_abs_diff_eq!(du, 0.9, epsilon = 1e-14);
    assert_abs_diff_eq!(d2u, 0.09, epsilon = 1e-14);
  }

  #[test]
  fn test_optimize_method_chain_rule_log_small_t() {
    let t = 1e-10;
    let (du, d2u) = chain_rule_log(t, 1e6, -1e12);
    assert_abs_diff_eq!(du, 1e-4, epsilon = 1e-14);
    assert_abs_diff_eq!(d2u, -1e-8 + 1e-4, epsilon = 1e-14);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::small(   0.01,  2, 10.0)]
  #[case::medium(  0.1,   3, 20.0)]
  #[case::large(   0.5,   1,  5.0)]
  #[case::no_indel(0.2,   0,  0.0)]
  #[trace]
  fn test_optimize_method_chain_rule_log_numerical_first_derivative(
    #[case] t: f64,
    #[case] k: usize,
    #[case] mu: f64,
  ) {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.5, 0.3, 0.1, 0.1]];
    let contribution = OptimizationContribution::Dense(
      optimize::dense::PartitionContribution::new(coefficients, gtr),
    );
    let contributions = vec![contribution];

    let u = t.ln();
    let metrics = evaluate_mixed(&contributions, t).expect("valid branch length");
    let indel = poisson_indel_log_lh(k, mu, t).expect("valid Poisson parameters");
    let dl_dt = metrics.derivative + indel.derivative;
    let d2l_dt2 = metrics.second_derivative + indel.second_derivative;
    let (dl_du_analytical, _) = chain_rule_log(t, dl_dt, d2l_dt2);

    let h = u.abs() * 1e-5;
    let eval_u = |uv: f64| {
      let tv = uv.exp();
      evaluate_mixed_log_lh_only(&contributions, tv).expect("valid branch length").value() + poisson_indel_log_lh(k, mu, tv).expect("valid Poisson parameters").log_lh.value()
    };
    let dl_du_numerical = (eval_u(u + h) - eval_u(u - h)) / (2.0 * h);

    assert_abs_diff_eq!(dl_du_analytical, dl_du_numerical, epsilon = 1e-4);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::small(   0.01,  2, 10.0)]
  #[case::medium(  0.1,   3, 20.0)]
  #[case::large(   0.5,   1,  5.0)]
  #[case::no_indel(0.2,   0,  0.0)]
  #[trace]
  fn test_optimize_method_chain_rule_log_numerical_second_derivative(
    #[case] t: f64,
    #[case] k: usize,
    #[case] mu: f64,
  ) {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.5, 0.3, 0.1, 0.1]];
    let contribution = OptimizationContribution::Dense(
      optimize::dense::PartitionContribution::new(coefficients, gtr),
    );
    let contributions = vec![contribution];

    let u = t.ln();
    let metrics = evaluate_mixed(&contributions, t).expect("valid branch length");
    let indel = poisson_indel_log_lh(k, mu, t).expect("valid Poisson parameters");
    let dl_dt = metrics.derivative + indel.derivative;
    let d2l_dt2 = metrics.second_derivative + indel.second_derivative;
    let (_, d2l_du2_analytical) = chain_rule_log(t, dl_dt, d2l_dt2);

    let h = u.abs() * 1e-4;
    let eval_u = |uv: f64| {
      let tv = uv.exp();
      evaluate_mixed_log_lh_only(&contributions, tv).expect("valid branch length").value() + poisson_indel_log_lh(k, mu, tv).expect("valid Poisson parameters").log_lh.value()
    };
    let d2l_du2_numerical = (eval_u(u + h) - 2.0 * eval_u(u) + eval_u(u - h)) / (h * h);

    assert_abs_diff_eq!(d2l_du2_analytical, d2l_du2_numerical, epsilon = 1e-2);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::small(   0.01,  2, 10.0)]
  #[case::medium(  0.1,   3, 20.0)]
  #[case::large(   0.5,   1,  5.0)]
  #[case::no_indel(0.2,   0,  0.0)]
  #[trace]
  fn test_optimize_method_chain_rule_numerical_first_derivative(
    #[case] t: f64,
    #[case] k: usize,
    #[case] mu: f64,
  ) {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.5, 0.3, 0.1, 0.1]];
    let contribution = OptimizationContribution::Dense(
      optimize::dense::PartitionContribution::new(coefficients, gtr),
    );
    let contributions = vec![contribution];

    let s = t.sqrt();
    let metrics = evaluate_mixed(&contributions, t).expect("valid branch length");
    let indel = poisson_indel_log_lh(k, mu, t).expect("valid Poisson parameters");
    let dl_dt = metrics.derivative + indel.derivative;
    let d2l_dt2 = metrics.second_derivative + indel.second_derivative;
    let (dl_ds_analytical, _) = chain_rule_sqrt(s, dl_dt, d2l_dt2);

    let h = s * 1e-5;
    let eval_s = |sv: f64| {
      let tv = sv * sv;
      evaluate_mixed_log_lh_only(&contributions, tv).expect("valid branch length").value() + poisson_indel_log_lh(k, mu, tv).expect("valid Poisson parameters").log_lh.value()
    };
    let dl_ds_numerical = (eval_s(s + h) - eval_s(s - h)) / (2.0 * h);

    assert_abs_diff_eq!(dl_ds_analytical, dl_ds_numerical, epsilon = 1e-4);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::small(   0.01,  2, 10.0)]
  #[case::medium(  0.1,   3, 20.0)]
  #[case::large(   0.5,   1,  5.0)]
  #[case::no_indel(0.2,   0,  0.0)]
  #[trace]
  fn test_optimize_method_chain_rule_numerical_second_derivative(
    #[case] t: f64,
    #[case] k: usize,
    #[case] mu: f64,
  ) {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.5, 0.3, 0.1, 0.1]];
    let contribution = OptimizationContribution::Dense(
      optimize::dense::PartitionContribution::new(coefficients, gtr),
    );
    let contributions = vec![contribution];

    let s = t.sqrt();
    let metrics = evaluate_mixed(&contributions, t).expect("valid branch length");
    let indel = poisson_indel_log_lh(k, mu, t).expect("valid Poisson parameters");
    let dl_dt = metrics.derivative + indel.derivative;
    let d2l_dt2 = metrics.second_derivative + indel.second_derivative;
    let (_, d2l_ds2_analytical) = chain_rule_sqrt(s, dl_dt, d2l_dt2);

    let h = s * 1e-4;
    let eval_s = |sv: f64| {
      let tv = sv * sv;
      evaluate_mixed_log_lh_only(&contributions, tv).expect("valid branch length").value() + poisson_indel_log_lh(k, mu, tv).expect("valid Poisson parameters").log_lh.value()
    };
    let d2l_ds2_numerical = (eval_s(s + h) - 2.0 * eval_s(s) + eval_s(s - h)) / (h * h);

    assert_abs_diff_eq!(d2l_ds2_analytical, d2l_ds2_numerical, epsilon = 1e-2);
  }
}
