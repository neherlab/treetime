#[cfg(test)]
mod tests {
  use crate::traits::ConvolveAlgo;
  use crate::{DirectConvolve, FftConvolve, NdarrayConvolve};
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

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
}
