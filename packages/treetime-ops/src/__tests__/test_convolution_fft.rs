#[cfg(test)]
mod tests {
  use crate::convolution::convolve_fft;
  use ndarray::Array1;
  use ndarray_conv::{ConvFFTExt, ConvMode, PaddingMode};

  #[test]
  fn test_convolution_fft_equals_fresh_plans_across_interleaved_sizes() {
    let lengths = helpers::interleaved_lengths();
    let operands: Vec<_> = lengths
      .iter()
      .enumerate()
      .map(|(seed, &(f_len, g_len))| helpers::operands(f_len, g_len, seed))
      .collect();

    let expected: Vec<Array1<f64>> = operands.iter().map(|(f, g)| helpers::fresh(0.1, f, g)).collect();
    let actual: Vec<Array1<f64>> = operands.iter().map(|(f, g)| convolve_fft(0.1, f, g).unwrap()).collect();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_fft_equals_fresh_plans_after_eviction() {
    let (f, g) = helpers::operands(1500, 700, 0);
    let expected = helpers::fresh(0.25, &f, &g);
    let first = convolve_fft(0.25, &f, &g).unwrap();
    for f_len in 100..200 {
      let (a, b) = helpers::operands(f_len * 7, f_len, f_len);
      convolve_fft(0.25, &a, &b).unwrap();
    }
    let after_eviction = convolve_fft(0.25, &f, &g).unwrap();
    assert_eq!((&expected, &expected), (&first, &after_eviction));
  }

  mod helpers {
    use super::*;

    pub(super) fn interleaved_lengths() -> Vec<(usize, usize)> {
      let sizes: Vec<(usize, usize)> = (100..1000).step_by(50).map(|g_len| (g_len * 3 + 7, g_len)).collect();
      sizes.iter().chain(sizes.iter().rev()).copied().collect()
    }

    pub(super) fn operands(f_len: usize, g_len: usize, seed: usize) -> (Array1<f64>, Array1<f64>) {
      (values(f_len, seed), values(g_len, seed + 1))
    }

    pub(super) fn fresh(dx: f64, f: &Array1<f64>, g: &Array1<f64>) -> Array1<f64> {
      &f.conv_fft(g, ConvMode::Full, PaddingMode::Zeros).unwrap() * dx
    }

    #[allow(clippy::as_conversions, reason = "small test indices convert exactly to f64")]
    fn values(len: usize, seed: usize) -> Array1<f64> {
      Array1::from_shape_fn(len, |i| ((i * 7 + seed * 13) as f64 * 0.37).sin().abs() + 0.01)
    }
  }
}
