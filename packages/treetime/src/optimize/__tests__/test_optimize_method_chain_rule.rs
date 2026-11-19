#[cfg(test)]
mod tests {

  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::indel::poisson_indel_log_lh;
  use crate::optimize::likelihood::{evaluate_partitions, evaluate_partitions_log_lh_only};
  use crate::optimize::method_newton::chain_rule_sqrt;
  use crate::partition::optimize;
  use crate::partition::optimize::contribution::OptimizationContribution;
  use ndarray::array;
  use proptest::prelude::*;
  use rstest::rstest;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[test]
  fn test_optimize_method_chain_rule_at_zero() {
    let (ds, d2s) = chain_rule_sqrt(0.0, 100.0, -500.0);
    pretty_assert_abs_diff_eq!(0.0, ds, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(200.0, d2s, epsilon = 1e-15);
  }

  #[test]
  fn test_optimize_method_chain_rule_analytical() {
    let (ds, d2s) = chain_rule_sqrt(0.3, 10.0, -100.0);
    pretty_assert_abs_diff_eq!(6.0, ds, epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(-16.0, d2s, epsilon = 1e-14);
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
      optimize::dense::PartitionContribution::new(coefficients, &gtr),
    );
    let contributions = vec![contribution];

    let s = t.sqrt();
    let metrics = evaluate_partitions(&contributions, t).expect("valid branch length");
    let indel = poisson_indel_log_lh(k, mu, t).expect("valid Poisson parameters");
    let dl_dt = metrics.derivative + indel.derivative;
    let d2l_dt2 = metrics.second_derivative + indel.second_derivative;
    let (dl_ds_analytical, _) = chain_rule_sqrt(s, dl_dt, d2l_dt2);

    let h = s * 1e-5;
    let eval_s = |sv: f64| {
      let tv = sv * sv;
      evaluate_partitions_log_lh_only(&contributions, tv).expect("valid branch length").value() + poisson_indel_log_lh(k, mu, tv).expect("valid Poisson parameters").log_lh.value()
    };
    let dl_ds_numerical = (eval_s(s + h) - eval_s(s - h)) / (2.0 * h);

    pretty_assert_abs_diff_eq!(dl_ds_analytical, dl_ds_numerical, epsilon = 1e-4);
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
      optimize::dense::PartitionContribution::new(coefficients, &gtr),
    );
    let contributions = vec![contribution];

    let s = t.sqrt();
    let metrics = evaluate_partitions(&contributions, t).expect("valid branch length");
    let indel = poisson_indel_log_lh(k, mu, t).expect("valid Poisson parameters");
    let dl_dt = metrics.derivative + indel.derivative;
    let d2l_dt2 = metrics.second_derivative + indel.second_derivative;
    let (_, d2l_ds2_analytical) = chain_rule_sqrt(s, dl_dt, d2l_dt2);

    let h = s * 1e-4;
    let eval_s = |sv: f64| {
      let tv = sv * sv;
      evaluate_partitions_log_lh_only(&contributions, tv).expect("valid branch length").value() + poisson_indel_log_lh(k, mu, tv).expect("valid Poisson parameters").log_lh.value()
    };
    let d2l_ds2_numerical = (eval_s(s + h) - 2.0 * eval_s(s) + eval_s(s - h)) / (h * h);

    pretty_assert_abs_diff_eq!(d2l_ds2_analytical, d2l_ds2_numerical, epsilon = 1e-2);
  }

  proptest! {
    #[test]
    fn test_prop_optimize_method_chain_rule_sqrt_formula(
      s in generators::gen_s(),
      dl_dt in generators::gen_dl_dt(),
      d2l_dt2 in generators::gen_d2l_dt2(),
    ) {
      let (dl_ds, d2l_ds2) = chain_rule_sqrt(s, dl_dt, d2l_dt2);
      let expected_dl_ds = 2.0 * s * dl_dt;
      let expected_d2l_ds2 = 4.0 * s * s * d2l_dt2 + 2.0 * dl_dt;
      let dl_tol = 1e-9 * expected_dl_ds.abs().max(1.0);
      let d2l_tol = 1e-9 * expected_d2l_ds2.abs().max(1.0);
      prop_assert!(
        (dl_ds - expected_dl_ds).abs() <= dl_tol,
        "dl_ds: got {dl_ds}, expected {expected_dl_ds}"
      );
      prop_assert!(
        (d2l_ds2 - expected_d2l_ds2).abs() <= d2l_tol,
        "d2l_ds2: got {d2l_ds2}, expected {expected_d2l_ds2}"
      );
    }

    #[test]
    fn test_prop_optimize_method_chain_rule_sqrt_linear(
      s in generators::gen_s(),
      dl_dt in generators::gen_dl_dt(),
      d2l_dt2 in generators::gen_d2l_dt2(),
      k in generators::gen_scalar(),
    ) {
      let (dl_ds, d2l_ds2) = chain_rule_sqrt(s, dl_dt, d2l_dt2);
      let (dl_ds_k, d2l_ds2_k) = chain_rule_sqrt(s, k * dl_dt, k * d2l_dt2);
      let dl_tol = 1e-9 * (k * dl_ds).abs().max(1.0);
      let d2l_tol = 1e-9 * (k * d2l_ds2).abs().max(1.0);
      prop_assert!(
        (dl_ds_k - k * dl_ds).abs() <= dl_tol,
        "linearity violated for dl_ds: {dl_ds_k} vs k*{dl_ds}"
      );
      prop_assert!(
        (d2l_ds2_k - k * d2l_ds2).abs() <= d2l_tol,
        "linearity violated for d2l_ds2: {d2l_ds2_k} vs k*{d2l_ds2}"
      );
    }
  }

  mod generators {
    use proptest::prelude::*;

    pub(super) fn gen_s() -> impl Strategy<Value = f64> {
      1e-6_f64..1e3_f64
    }

    pub(super) fn gen_dl_dt() -> impl Strategy<Value = f64> {
      -1e6_f64..1e6_f64
    }

    pub(super) fn gen_d2l_dt2() -> impl Strategy<Value = f64> {
      -1e8_f64..1e8_f64
    }

    pub(super) fn gen_scalar() -> impl Strategy<Value = f64> {
      prop_oneof![-1e3_f64..-1e-3_f64, 1e-3_f64..1e3_f64]
    }
  }
}
