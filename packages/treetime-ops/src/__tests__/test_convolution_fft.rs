#[cfg(test)]
mod tests {
  use crate::convolution::{FftConvolver, FftPlans, convolve_fft, smooth_even_len};
  use helpers::Shape;
  use itertools::Itertools;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
  use rstest::rstest;
  use std::iter::repeat_with;
  use std::sync::Arc;
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::gaussian_smallest(                Shape::Gaussian, (   2,    2))]
  #[case::noise_smallest(                   Shape::Noise,    (   2,    2))]
  #[case::gaussian_span_at_smallest_length( Shape::Gaussian, (   3,    2))]
  #[case::gaussian_odd_span(                Shape::Gaussian, (1600,   40))]
  #[case::noise_odd_span(                   Shape::Noise,    (1600,   40))]
  #[case::gaussian_even_span(               Shape::Gaussian, (  37, 1200))]
  #[case::noise_even_span(                  Shape::Noise,    (  37, 1200))]
  #[case::gaussian_long_second_operand(     Shape::Gaussian, (   5, 4000))]
  #[case::noise_long_second_operand(        Shape::Noise,    (   5, 4000))]
  #[case::gaussian_span_at_smooth_length(   Shape::Gaussian, (1689,   40))]
  #[case::noise_span_at_smooth_length(      Shape::Noise,    (1689,   40))]
  #[case::gaussian_span_above_smooth_length(Shape::Gaussian, (1690,   40))]
  #[case::noise_span_above_smooth_length(   Shape::Noise,    (1690,   40))]
  #[case::gaussian_span_at_power_of_two(    Shape::Gaussian, (2000,   49))]
  #[case::noise_span_above_power_of_two(    Shape::Noise,    (2000,   50))]
  #[case::gaussian_largest_odd_span(        Shape::Gaussian, (4000, 4000))]
  #[case::noise_largest_odd_span(           Shape::Noise,    (4000, 4000))]
  #[case::gaussian_largest_even_span(       Shape::Gaussian, (4000, 3999))]
  #[case::noise_largest_even_span(          Shape::Noise,    (4000, 3999))]
  #[trace]
  fn test_convolution_fft_equals_the_direct_sum_within_rounding_of_its_peak(
    #[case] shape: Shape,
    #[case] (f_len, g_len): (usize, usize),
  ) {
    let (f, g) = (shape.values(f_len, 0), shape.values(g_len, 1));
    let direct = helpers::direct(0.1, &f, &g);
    let actual = convolve_fft(0.1, &f, &g).unwrap();
    pretty_assert_abs_diff_eq!(
      helpers::over_peak_of(&direct, &direct),
      helpers::over_peak_of(&direct, &actual),
      epsilon = 1e-14
    );
  }

  #[test]
  fn test_convolution_fft_equals_the_direct_sum_for_operand_pairs_of_one_transform_length() {
    let operands = helpers::operands_of(&[(100, 900), (950, 50), (300, 700), (700, 300), (100, 900), (950, 50)]);

    let plans = FftPlans::new();
    let actual = operands
      .iter()
      .map(|(f, g)| {
        FftConvolver::new(f.len(), g.len(), &plans)
          .unwrap()
          .convolve(0.1, f, g)
          .unwrap()
      })
      .collect_vec();

    helpers::assert_equal_to_direct_sums(0.1, &operands, &actual);
  }

  #[test]
  fn test_convolution_fft_convolver_reused_on_new_operands_equals_the_direct_sum() {
    let plans = FftPlans::new();
    let mut reused = FftConvolver::new(1500, 700, &plans).unwrap();
    let rounds = (0..4)
      .map(|round| {
        let (f, g) = helpers::operands(1500, 700, 10 * round);
        let (a, b) = helpers::operands(300 + 100 * round, 40 + 10 * round, 10 * round + 5);
        (f, g, a, b)
      })
      .collect_vec();

    let (operands, actual): (Vec<_>, Vec<_>) = rounds
      .iter()
      .flat_map(|(f, g, a, b)| {
        let reused_result = reused.convolve(0.25, f, g).unwrap();
        let other_result = FftConvolver::new(a.len(), b.len(), &plans)
          .unwrap()
          .convolve(0.25, a, b)
          .unwrap();
        [
          ((f.clone(), g.clone()), reused_result),
          ((a.clone(), b.clone()), other_result),
        ]
      })
      .unzip();

    helpers::assert_equal_to_direct_sums(0.25, &operands, &actual);
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
  fn test_convolution_fft_shared_plans_in_parallel_equal_the_direct_sum() {
    let lengths = helpers::interleaved_lengths()
      .into_iter()
      .chain([(1600, 40), (1600, 40), (37, 1200), (2, 2), (1689, 40), (5000, 3000)])
      .collect_vec();
    let operands = helpers::operands_of(&lengths);

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

    helpers::assert_equal_to_direct_sums(0.1, &operands, &actual);
  }

  #[test]
  fn test_convolution_fft_plans_of_one_length_are_shared_by_every_lookup() {
    let plans = FftPlans::new();
    let concurrent = helpers::pool().install(|| (0..64).into_par_iter().map(|_| plans.get(1728)).collect::<Vec<_>>());
    let sequential = repeat_with(|| plans.get(1728)).take(2).collect_vec();

    let first = &concurrent[0];
    let shares_first = concurrent
      .iter()
      .chain(&sequential)
      .map(|pair| Arc::ptr_eq(&pair.forward, &first.forward) && Arc::ptr_eq(&pair.inverse, &first.inverse))
      .collect_vec();
    assert_eq!(vec![true; 66], shares_first);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::below_four(            3,    4)]
  #[case::below_eight(           7,    8)]
  #[case::below_twelve(         11,   12)]
  #[case::skips_odd_smooth(     13,   16)]
  #[case::factors_two_and_three(1639, 1728)]
  #[case::power_of_two(       2048, 2048)]
  #[case::factors_of_all_three(2049, 2160)]
  #[trace]
  fn test_convolution_fft_transform_length_is_the_smallest_even_five_smooth_number_not_below_the_span(
    #[case] span: usize,
    #[case] expected: usize,
  ) {
    assert_eq!(expected, smooth_even_len(span));
  }

  mod helpers {
    use crate::DirectConvolve;
    use crate::traits::ConvolveAlgo;
    use itertools::{Itertools, izip};
    use ndarray::Array1;
    use rayon::{ThreadPool, ThreadPoolBuilder};
    use treetime_utils::array::ndarray::max_or;
    use treetime_utils::cast::ToF64Exact;
    use treetime_utils::pretty_assert_abs_diff_eq;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Shape {
      Gaussian,
      Noise,
    }

    impl Shape {
      pub(super) fn values(self, len: usize, seed: usize) -> Array1<f64> {
        let len_f64 = len.to_f64_exact();
        match self {
          Self::Gaussian => {
            let (center, width) = (len_f64 * (0.4 + 0.1 * seed.to_f64_exact()), len_f64 / 8.0 + 1.0);
            Array1::from_shape_fn(len, |i| (-0.5 * ((i.to_f64_exact() - center) / width).powi(2)).exp())
          },
          Self::Noise => Array1::from_shape_fn(len, |i| {
            ((i + 7919 * seed).to_f64_exact() * 0.618_033_988_749_894_9).fract() + 0.01
          }),
        }
      }
    }

    pub(super) fn direct(dx: f64, f: &Array1<f64>, g: &Array1<f64>) -> Array1<f64> {
      DirectConvolve.convolve(dx, f, g).unwrap()
    }

    pub(super) fn over_peak_of(reference: &Array1<f64>, values: &Array1<f64>) -> Array1<f64> {
      values / max_or(reference, f64::NAN)
    }

    pub(super) fn assert_equal_to_direct_sums(
      dx: f64,
      operands: &[(Array1<f64>, Array1<f64>)],
      actual: &[Array1<f64>],
    ) {
      let direct_sums = operands.iter().map(|(f, g)| direct(dx, f, g)).collect_vec();
      let expected = direct_sums.iter().map(|d| over_peak_of(d, d)).collect_vec();
      let actual = izip!(&direct_sums, actual)
        .map(|(d, a)| over_peak_of(d, a))
        .collect_vec();
      pretty_assert_abs_diff_eq!(expected[..], actual[..], epsilon = 1e-14);
    }

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

    pub(super) fn pool() -> ThreadPool {
      ThreadPoolBuilder::new().num_threads(8).build().unwrap()
    }

    fn values(len: usize, seed: usize) -> Array1<f64> {
      Array1::from_shape_fn(len, |i| ((i * 7 + seed * 13).to_f64_exact() * 0.37).sin().abs() + 0.01)
    }
  }
}
