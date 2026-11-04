#[cfg(test)]
mod tests {
  use crate::traits::ConvolveAlgo;
  use crate::{
    DirectConvolve, FftConvolve, FftConvolver, FftPlans, NdarrayConvolve, convolve_direct_at,
    convolve_direct_neg_log_at,
  };
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  #[rstest]
  #[case::direct(&DirectConvolve)]
  #[case::ndarray_conv(&NdarrayConvolve)]
  #[case::fft(&FftConvolve)]
  #[trace]
  fn test_convolution_algorithms_agree_with_the_trapezoid_rule_at_hard_ends(#[case] algo: &dyn ConvolveAlgo) {
    let (f, g) = (array![1.0, 2.0, 1.0], array![1.0, 2.0]);
    let trapezoid_by_hand = array![0.25, 2.0, 2.5, 0.5];
    let actual = helpers::scaled_to_unit_peak(&algo.convolve(1.0, &f, &g).unwrap());
    pretty_assert_abs_diff_eq!(
      helpers::scaled_to_unit_peak(&trapezoid_by_hand),
      actual,
      epsilon = 1e-15
    );

    let (f, g) = (helpers::exponential(1.5, 81, 0.05), helpers::exponential(4.0, 41, 0.05));
    let direct = helpers::scaled_to_unit_peak(&DirectConvolve.convolve(0.05, &f, &g).unwrap());
    let actual = helpers::scaled_to_unit_peak(&algo.convolve(0.05, &f, &g).unwrap());
    pretty_assert_abs_diff_eq!(direct, actual, epsilon = 1e-12);
  }

  #[rstest]
  #[case::direct(&DirectConvolve)]
  #[case::ndarray_conv(&NdarrayConvolve)]
  #[case::fft(&FftConvolve)]
  #[trace]
  fn test_convolution_inner_outputs_are_trapezoid_sums_and_extreme_outputs_are_quarter_end_products(
    #[case] algo: &dyn ConvolveAlgo,
  ) {
    let (f, g) = (array![2.0, 1.0, 3.0], array![5.0, 1.0, 1.0, 7.0]);
    let trapezoid_by_hand = array![1.25, 1.75, 4.75, 4.75, 2.5, 2.625];
    let actual = algo.convolve(0.5, &f, &g).unwrap();
    pretty_assert_abs_diff_eq!(trapezoid_by_hand, actual, epsilon = 1e-14);
  }

  #[test]
  fn test_convolve_direct_at_gives_each_trapezoid_output_alone() {
    let (f, g) = (array![2.0, 1.0, 3.0], array![5.0, 1.0, 1.0, 7.0]);
    let trapezoid_by_hand = array![1.25, 1.75, 4.75, 4.75, 2.5, 2.625];
    let actual: Array1<f64> = (0..trapezoid_by_hand.len())
      .map(|n| convolve_direct_at(0.5, &f, &g, n).unwrap())
      .collect();
    pretty_assert_abs_diff_eq!(trapezoid_by_hand, actual, epsilon = 1e-14);
    assert_error!(
      convolve_direct_at(0.5, &f, &g, 6),
      "Convolution sample 6 is outside the full convolution of 6 samples"
    );
  }

  #[test]
  fn test_convolution_overlap_sums_are_the_plain_sums_of_products_without_trapezoid_weights() {
    let (f, g) = (array![2.0, 1.0, 3.0], array![5.0, 1.0, 1.0, 7.0]);
    let plain_sums_by_hand = array![10.0, 7.0, 18.0, 18.0, 10.0, 21.0];
    let actual = FftConvolver::new(3, 4, &FftPlans::new())
      .unwrap()
      .overlap_sums(&f, &g)
      .unwrap();
    pretty_assert_abs_diff_eq!(plain_sums_by_hand, actual, epsilon = 1e-13);
  }

  #[test]
  fn test_convolve_direct_neg_log_at_is_the_neg_log_of_the_trapezoid_sum() {
    let (f, g) = (array![2.0, 1.0, 3.0], array![5.0, 1.0, 1.0, 7.0]);
    let (f_neg_log, g_neg_log) = (f.mapv(|v: f64| -v.ln()), g.mapv(|v: f64| -v.ln()));
    let expected: Array1<f64> = (0..6)
      .map(|n| -convolve_direct_at(0.5, &f, &g, n).unwrap().ln())
      .collect();
    let actual: Array1<f64> = (0..6)
      .map(|n| convolve_direct_neg_log_at(0.5, &f_neg_log, &g_neg_log, n, 0..=n).unwrap())
      .collect();
    pretty_assert_abs_diff_eq!(expected, actual, epsilon = 1e-14);
    assert_error!(
      convolve_direct_neg_log_at(0.5, &f_neg_log, &g_neg_log, 6, 0..=6),
      "Convolution sample 6 is outside the full convolution of 6 samples"
    );
  }

  #[test]
  fn test_convolve_direct_neg_log_at_rejects_a_nan_operand_value_in_its_pairs() {
    let (f_neg_log, g_neg_log) = (array![1.0, f64::NAN, 2.0], array![1.0, 3.0]);
    assert_error!(
      convolve_direct_neg_log_at(0.5, &f_neg_log, &g_neg_log, 2, 0..=2),
      "Convolution sample 2 has a NaN operand value in its pair (1, 1)"
    );
  }

  #[test]
  fn test_convolve_direct_neg_log_at_keeps_depths_below_the_smallest_double_and_exact_zeros() {
    let f = array![1000.0, 1000.0, f64::INFINITY];
    let g = array![2000.0, 2000.0, f64::INFINITY];
    let at = |n: usize| convolve_direct_neg_log_at(0.5, &f, &g, n, 0..=n).unwrap();
    pretty_assert_ulps_eq!(3000.0 - (0.5_f64).ln(), at(1), max_ulps = 4);
    assert_eq!(f64::INFINITY, at(4));
  }

  #[test]
  fn test_convolve_direct_neg_log_at_sums_only_the_pairs_of_its_index_range_with_the_weights_of_the_full_sum() {
    let (f, g) = (array![2.0, 1.0, 3.0], array![5.0, 1.0, 1.0, 7.0]);
    let (f_neg_log, g_neg_log) = (f.mapv(|v: f64| -v.ln()), g.mapv(|v: f64| -v.ln()));
    let middle_pair_of_sample_2 = f[1] * g[1] * 0.5;
    let actual = convolve_direct_neg_log_at(0.5, &f_neg_log, &g_neg_log, 2, 1..=1).unwrap();
    pretty_assert_abs_diff_eq!(-middle_pair_of_sample_2.ln(), actual, epsilon = 1e-15);
    let end_pairs_of_sample_2 = f64::midpoint(f[0] * g[2], f[2] * g[0]) * 0.5;
    let actual = convolve_direct_neg_log_at(0.5, &f_neg_log, &g_neg_log, 2, 0..=2).unwrap();
    pretty_assert_abs_diff_eq!(
      -(end_pairs_of_sample_2 + middle_pair_of_sample_2).ln(),
      actual,
      epsilon = 1e-15
    );
    assert_eq!(
      f64::INFINITY,
      convolve_direct_neg_log_at(0.5, &f_neg_log, &g_neg_log, 2, 3..=5).unwrap()
    );
  }

  #[test]
  fn test_convolution_direct_delta_shifts_the_operand_with_half_weight_at_its_grid_ends() {
    let dx = 0.1;
    let delta = array![0.0, 0.0, 1.0 / dx, 0.0, 0.0];
    let f = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let shifted_by_hand = array![0.0, 0.0, 0.5, 2.0, 3.0, 2.0, 0.5, 0.0, 0.0];
    let actual = DirectConvolve.convolve(dx, &delta, &f).unwrap();
    pretty_assert_abs_diff_eq!(shifted_by_hand, actual, epsilon = 1e-15);
  }

  #[test]
  fn test_convolution_is_symmetric_in_its_operands() {
    let f = array![1.0, 2.0, 1.0];
    let g = array![1.0, 1.0, 1.0];
    let forward = NdarrayConvolve.convolve(0.1, &f, &g).unwrap();
    let backward = NdarrayConvolve.convolve(0.1, &g, &f).unwrap();
    pretty_assert_ulps_eq!(forward, backward, max_ulps = 4);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::direct_one_sample(      &DirectConvolve,  (array![2.0],      array![1.0, 2.0]), "got 1 and 2")]
  #[case::direct_empty(           &DirectConvolve,  (array![1.0, 2.0], Array1::zeros(0)), "got 2 and 0")]
  #[case::ndarray_conv_one_sample(&NdarrayConvolve, (array![2.0],      array![1.0, 2.0]), "got 1 and 2")]
  #[case::ndarray_conv_empty(     &NdarrayConvolve, (array![1.0, 2.0], Array1::zeros(0)), "got 2 and 0")]
  #[case::fft_one_sample(         &FftConvolve,     (array![2.0],      array![1.0, 2.0]), "got 1 and 2")]
  #[case::fft_empty(              &FftConvolve,     (array![1.0, 2.0], Array1::zeros(0)), "got 2 and 0")]
  #[trace]
  fn test_convolution_rejects_operands_with_fewer_than_two_samples(
    #[case] algo: &dyn ConvolveAlgo,
    #[case] (f, g): (Array1<f64>, Array1<f64>),
    #[case] lengths: &str,
  ) {
    assert_error!(
      algo.convolve(1.0, &f, &g),
      format!("Trapezoid convolution needs at least 2 samples per operand, {lengths}")
    );
  }

  mod helpers {
    use ndarray::Array1;
    use treetime_utils::array::ndarray::max_or;
    use treetime_utils::cast::ToF64Exact;

    pub(super) fn exponential(rate: f64, len: usize, dx: f64) -> Array1<f64> {
      Array1::from_shape_fn(len, |i| (-rate * i.to_f64_exact() * dx).exp())
    }

    pub(super) fn scaled_to_unit_peak(values: &Array1<f64>) -> Array1<f64> {
      values / max_or(values, f64::NAN)
    }
  }
}
