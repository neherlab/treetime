#[cfg(test)]
mod tests {
  use crate::multiplication::{multiply_many, multiply_many_lazy_normalize};
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::iter::repeat_n;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::{pretty_assert_neg_inf, pretty_assert_ulps_eq};

  #[test]
  fn test_multiply_many_lazy_normalize_of_no_factors_is_empty_and_impossible() {
    let result = multiply_many_lazy_normalize(&[]);
    assert_eq!(Array1::<f64>::zeros(0), result.normalized);
    pretty_assert_neg_inf!(result.log_scale);
  }

  #[test]
  fn test_multiply_many_lazy_normalize_scales_single_factor_to_unit_maximum() {
    let a = array![1.0, 4.0, 2.0];
    let result = multiply_many_lazy_normalize(&[&a]);

    pretty_assert_ulps_eq!(array![0.25, 1.0, 0.5], result.normalized, max_ulps = 4);
    pretty_assert_ulps_eq!(4.0_f64.ln(), result.log_scale, max_ulps = 4);
  }

  #[test]
  fn test_multiply_many_lazy_normalize_keeps_product_as_normalized_times_scale() {
    let a = array![1.0, 2.0, 1.0];
    let result = multiply_many_lazy_normalize(&[&a, &a]);

    pretty_assert_ulps_eq!(array![0.25, 1.0, 0.25], result.normalized, max_ulps = 4);
    pretty_assert_ulps_eq!(4.0_f64.ln(), result.log_scale, max_ulps = 4);
  }

  #[test]
  fn test_multiply_many_lazy_normalize_moves_underflowing_product_into_log_scale() {
    let n = 100;
    let small_val = 0.0001;
    let a = array![small_val];

    let refs: Vec<&Array1<f64>> = repeat_n(&a, n).collect();
    let result = multiply_many_lazy_normalize(&refs);

    let expected_log_scale = n.to_f64_exact() * small_val.ln();
    pretty_assert_ulps_eq!(array![1.0], result.normalized, max_ulps = 4);
    pretty_assert_ulps_eq!(expected_log_scale, result.log_scale, max_ulps = 4);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::single_zero_factor(vec![array![0.0, 0.0, 0.0]])]
  #[case::zero_first_factor( vec![array![0.0, 0.0, 0.0], array![1.0, 1.0, 1.0]])]
  #[case::disjoint_factors(  vec![array![1.0, 0.0, 0.0], array![0.0, 1.0, 1.0]])]
  #[trace]
  fn test_multiply_many_lazy_normalize_of_factors_without_common_support_is_impossible(
    #[case] factors: Vec<Array1<f64>>,
  ) {
    let factors: Vec<&Array1<f64>> = factors.iter().collect();
    let result = multiply_many_lazy_normalize(&factors);

    pretty_assert_neg_inf!(result.log_scale);
    assert_eq!(Array1::<f64>::zeros(3), result.normalized);
  }

  #[test]
  fn test_multiply_many_lazy_normalize_matches_aggressive_normalization() {
    let a = array![1.0, 2.0, 1.0];
    let b = array![0.5, 1.0, 0.5];

    let lazy = multiply_many_lazy_normalize(&[&a, &b]);
    let agg = multiply_many(&[&a, &b]);

    pretty_assert_ulps_eq!(agg.normalized, lazy.normalized, max_ulps = 4);
    pretty_assert_ulps_eq!(agg.log_scale, lazy.log_scale, max_ulps = 4);
  }

  #[test]
  fn test_multiply_many_lazy_normalize_of_unit_maximum_factors_has_zero_scale() {
    let a = array![0.5, 1.0, 0.5];

    let result = multiply_many_lazy_normalize(&[&a, &a, &a]);

    pretty_assert_ulps_eq!(array![0.125, 1.0, 0.125], result.normalized, max_ulps = 4);
    pretty_assert_ulps_eq!(0.0, result.log_scale, max_ulps = 4);
  }
}
