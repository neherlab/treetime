#[cfg(test)]
mod tests {
  use crate::convolution::{FftConvolver, FftPlans, convolve_fft};
  use itertools::Itertools;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
  use rstest::rstest;
  use std::sync::Arc;
  use treetime_utils::assert_error;

  #[test]
  fn test_convolution_fft_equals_ndarray_conv_across_interleaved_sizes() {
    let operands = helpers::operands_of(&helpers::interleaved_lengths());

    let expected = operands.iter().map(|(f, g)| helpers::fresh(0.1, f, g)).collect_vec();
    let actual = operands
      .iter()
      .map(|(f, g)| convolve_fft(0.1, f, g).unwrap())
      .collect_vec();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_fft_equals_ndarray_conv_for_equal_output_lengths_of_different_padding() {
    let operands = helpers::operands_of(&[(100, 900), (950, 50), (300, 700), (700, 300), (100, 900), (950, 50)]);

    let expected = operands.iter().map(|(f, g)| helpers::fresh(0.1, f, g)).collect_vec();
    let actual = operands
      .iter()
      .map(|(f, g)| convolve_fft(0.1, f, g).unwrap())
      .collect_vec();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_fft_equals_ndarray_conv_for_an_odd_transform_length() {
    let (f, g) = helpers::operands(1600, 40, 3);
    let expected = helpers::fresh(0.2, &f, &g);
    let actual = convolve_fft(0.2, &f, &g).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_fft_convolver_reused_on_new_operands_equals_ndarray_conv() {
    let plans = FftPlans::new();
    let mut reused = FftConvolver::new(1500, 700, &plans).unwrap();
    let rounds = (0..4)
      .map(|round| {
        let (f, g) = helpers::operands(1500, 700, 10 * round);
        let (a, b) = helpers::operands(300 + 100 * round, 40 + 10 * round, 10 * round + 5);
        (f, g, a, b)
      })
      .collect_vec();

    let expected = rounds
      .iter()
      .map(|(f, g, a, b)| (helpers::fresh(0.25, f, g), helpers::fresh(0.25, a, b)))
      .collect_vec();
    let actual = rounds
      .iter()
      .map(|(f, g, a, b)| {
        let reused_result = reused.convolve(0.25, f, g).unwrap();
        let other_result = FftConvolver::new(a.len(), b.len(), &plans)
          .unwrap()
          .convolve(0.25, a, b)
          .unwrap();
        (reused_result, other_result)
      })
      .collect_vec();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_fft_convolver_rejects_operands_of_other_lengths() {
    let plans = FftPlans::new();
    let mut convolver = FftConvolver::new(5, 3, &plans).unwrap();
    let (f, g) = (array![1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0]);
    assert_error!(
      convolver.convolve(0.1, &f, &g),
      "FFT convolver for operands of 5 and 3 samples got operands of 4 and 3 samples"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::empty_first( (0, 3))]
  #[case::single_first((1, 3))]
  #[case::single_second((3, 1))]
  #[case::empty_second((3, 0))]
  #[trace]
  fn test_convolution_fft_convolver_rejects_operands_below_two_samples(#[case] (f_len, g_len): (usize, usize)) {
    let plans = FftPlans::new();
    assert_error!(
      FftConvolver::new(f_len, g_len, &plans),
      format!("Trapezoid convolution needs at least 2 samples per operand, got {f_len} and {g_len}")
    );
  }

  #[test]
  fn test_convolution_fft_shared_plans_in_parallel_equal_ndarray_conv() {
    let lengths = helpers::interleaved_lengths()
      .into_iter()
      .chain([(1600, 40), (1600, 40), (37, 1200), (2, 2), (5000, 3000)])
      .collect_vec();
    let operands = helpers::operands_of(&lengths);
    let expected = operands.iter().map(|(f, g)| helpers::fresh(0.1, f, g)).collect_vec();

    let plans = FftPlans::new();
    let actual = helpers::pool().install(|| {
      operands
        .par_iter()
        .map(|(f, g)| {
          FftConvolver::new(f.len(), g.len(), &plans)
            .unwrap()
            .convolve(0.1, f, g)
            .unwrap()
        })
        .collect::<Vec<_>>()
    });

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_fft_plans_of_one_length_are_shared_by_every_lookup() {
    let plans = FftPlans::new();
    let concurrent = helpers::pool().install(|| (0..64).into_par_iter().map(|_| plans.get(1705)).collect::<Vec<_>>());
    let sequential = (0..2).map(|_| plans.get(1705)).collect_vec();

    let first = &concurrent[0];
    let shares_first = concurrent
      .iter()
      .chain(&sequential)
      .map(|pair| Arc::ptr_eq(&pair.forward, &first.forward) && Arc::ptr_eq(&pair.inverse, &first.inverse))
      .collect_vec();
    assert_eq!(vec![true; 66], shares_first);
  }

  mod helpers {
    use crate::convolution::trapezoid_from_overlap_sums;
    use itertools::Itertools;
    use ndarray::Array1;
    use ndarray_conv::{ConvFFTExt, ConvMode, PaddingMode};
    use rayon::{ThreadPool, ThreadPoolBuilder};
    use treetime_utils::cast::ToF64Exact;

    pub(super) fn interleaved_lengths() -> Vec<(usize, usize)> {
      let sizes = (100..1000)
        .step_by(50)
        .map(|g_len| (g_len * 3 + 7, g_len))
        .collect_vec();
      sizes.iter().chain(sizes.iter().rev()).copied().collect()
    }

    pub(super) fn operands_of(lengths: &[(usize, usize)]) -> Vec<(Array1<f64>, Array1<f64>)> {
      lengths
        .iter()
        .enumerate()
        .map(|(seed, &(f_len, g_len))| operands(f_len, g_len, seed))
        .collect()
    }

    pub(super) fn operands(f_len: usize, g_len: usize, seed: usize) -> (Array1<f64>, Array1<f64>) {
      (values(f_len, seed), values(g_len, seed + 1))
    }

    pub(super) fn fresh(dx: f64, f: &Array1<f64>, g: &Array1<f64>) -> Array1<f64> {
      trapezoid_from_overlap_sums(f.conv_fft(g, ConvMode::Full, PaddingMode::Zeros).unwrap(), dx, f, g)
    }

    pub(super) fn pool() -> ThreadPool {
      ThreadPoolBuilder::new().num_threads(8).build().unwrap()
    }

    fn values(len: usize, seed: usize) -> Array1<f64> {
      Array1::from_shape_fn(len, |i| ((i * 7 + seed * 13).to_f64_exact() * 0.37).sin().abs() + 0.01)
    }
  }
}
