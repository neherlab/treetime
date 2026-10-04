use crate::traits::ConvolveAlgo;
use eyre::Report;
use ndarray::Array1;
use ndarray_conv::{ConvExt, ConvFFTExt, ConvMode, FftProcessor, PaddingMode};
use std::cell::RefCell;
use std::cmp::Ordering;

const FFT_PROCESSORS_PER_THREAD: usize = 8;
const MAX_CACHED_FFT_LEN: usize = 1 << 14;

thread_local! {
  static FFT_PROCESSORS_BY_FFT_LEN: RefCell<Vec<(usize, FftProcessor<f64>)>> = const { RefCell::new(Vec::new()) };
}

pub struct RiemannConvolve;

impl ConvolveAlgo for RiemannConvolve {
  fn name(&self) -> &'static str {
    "riemann"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve_riemann(dx, f_values, g_values)
  }
}

fn convolve_riemann(dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
  let mut result = Array1::zeros(f_values.len() + g_values.len() - 1);

  for (i, &f_val) in f_values.iter().enumerate() {
    for (j, &g_val) in g_values.iter().enumerate() {
      result[i + j] += f_val * g_val * dx;
    }
  }

  Ok(result)
}

pub struct NdarrayConvolve;

impl ConvolveAlgo for NdarrayConvolve {
  fn name(&self) -> &'static str {
    "ndarray-conv"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve(dx, f_values, g_values)
  }
}

fn convolve(dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
  let discrete_conv = f_values.conv(g_values, ConvMode::Full, PaddingMode::Zeros)?;
  let continuous_conv = &discrete_conv * dx;
  Ok(continuous_conv)
}

pub struct FftConvolve;

impl ConvolveAlgo for FftConvolve {
  fn name(&self) -> &'static str {
    "ndarray-conv-fft"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve_fft(dx, f_values, g_values)
  }
}

pub fn convolve_fft(dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
  let padded_len = f_values.len() + 2 * g_values.len().saturating_sub(1);
  let fft_len = fft_len(padded_len.max(g_values.len()));
  let discrete_conv = if fft_len > MAX_CACHED_FFT_LEN {
    f_values.conv_fft(g_values, ConvMode::Full, PaddingMode::Zeros)?
  } else {
    FFT_PROCESSORS_BY_FFT_LEN.with_borrow_mut(|processors| {
      let processor = processor_for_fft_len(processors, fft_len);
      f_values.conv_fft_with_processor(g_values, ConvMode::Full, PaddingMode::Zeros, processor)
    })?
  };
  let continuous_conv = &discrete_conv * dx;
  Ok(continuous_conv)
}

#[expect(
  clippy::integer_division,
  reason = "the transform length of ndarray-conv steps down from a power of two by integer ratios"
)]
fn fft_len(padded_len: usize) -> usize {
  let mut len = padded_len.next_power_of_two();
  for (numerator, denominator) in [(3, 4), (5, 6)] {
    loop {
      let smaller = len / denominator * numerator;
      match smaller.cmp(&padded_len) {
        Ordering::Less => break,
        Ordering::Equal => return padded_len,
        Ordering::Greater => len = smaller,
      }
    }
  }
  len
}

fn processor_for_fft_len(processors: &mut Vec<(usize, FftProcessor<f64>)>, fft_len: usize) -> &mut FftProcessor<f64> {
  let entry = if let Some(index) = processors.iter().position(|(len, _)| *len == fft_len) {
    processors.remove(index)
  } else {
    if processors.len() >= FFT_PROCESSORS_PER_THREAD {
      processors.remove(0);
    }
    (fft_len, FftProcessor::default())
  };
  processors.push(entry);
  let index = processors.len() - 1;
  &mut processors[index].1
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_convolve_riemann_delta() {
    let dx = 0.1;
    let delta = array![0.0, 0.0, 1.0 / dx, 0.0, 0.0];
    let f = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let result = convolve_riemann(dx, &delta, &f).unwrap();
    assert_eq!(9, result.len());
    assert_ulps_eq!(result[2], f[0], max_ulps = 4);
    assert_ulps_eq!(result[4], f[2], max_ulps = 4);
  }

  #[test]
  fn test_convolve_symmetric() {
    let dx = 0.1;
    let f = array![1.0, 2.0, 1.0];
    let g = array![1.0, 1.0, 1.0];
    let result1 = convolve(dx, &f, &g).unwrap();
    let result2 = convolve(dx, &g, &f).unwrap();
    for i in 0..result1.len() {
      assert_ulps_eq!(result1[i], result2[i], max_ulps = 4);
    }
  }

  #[test]
  fn test_convolve_fft_matches_direct() {
    let dx = 0.1;
    let f = array![1.0, 2.0, 3.0, 2.0, 1.0];
    let g = array![1.0, 1.0, 1.0];
    let result_direct = convolve(dx, &f, &g).unwrap();
    let result_fft = convolve_fft(dx, &f, &g).unwrap();
    for i in 0..result_direct.len() {
      assert_ulps_eq!(result_direct[i], result_fft[i], max_ulps = 100);
    }
  }
}
