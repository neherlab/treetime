use ndarray::Array1;
use statrs::function::erf::erfc;

pub fn gaussian_exponential_convolution_grid(a: f64, grid: &Array1<f64>) -> Array1<f64> {
  grid.mapv(|x| gaussian_exponential_convolution(a, x))
}

pub fn gaussian_exponential_convolution(a: f64, x: f64) -> f64 {
  0.5 * a * (-x * a + 0.5 * a.powi(2)).exp() * erfc((a - x) / 2_f64.sqrt())
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_gaussian_exponential_convolution_positive_x() {
    let result = gaussian_exponential_convolution(0.5, 1.0);
    assert_ulps_eq!(result, 0.23761736816002352, max_ulps = 4);
  }

  #[test]
  fn test_gaussian_exponential_convolution_grid() {
    let grid = array![-2.0, 0.0, 2.0, 5.0];
    let a = 0.5;
    let result = gaussian_exponential_convolution_grid(a, &grid);

    assert!(result[0] >= 0.0);
    assert!(result[1] >= 0.0);
    assert!(result[2] >= 0.0);
    assert!(result[3] >= 0.0);
  }

  #[test]
  fn test_gaussian_exponential_convolution_at_peak() {
    let a = 0.5;
    let result = gaussian_exponential_convolution(a, a);
    assert!(result > 0.0);
  }
}
