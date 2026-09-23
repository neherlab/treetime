use crate::traits::ConvolveAlgo;
use eyre::Report;
use ndarray::Array1;
use ndarray_conv::{ConvExt, ConvFFTExt, ConvMode, PaddingMode};

pub struct RiemannConvolve;

impl ConvolveAlgo for RiemannConvolve {
  fn name(&self) -> &'static str {
    "riemann"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve_riemann(dx, f_values, g_values)
  }
}

pub fn convolve_riemann(dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
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

pub fn convolve(dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
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
  let discrete_conv = f_values.conv_fft(g_values, ConvMode::Full, PaddingMode::Zeros)?;
  let continuous_conv = &discrete_conv * dx;
  Ok(continuous_conv)
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
