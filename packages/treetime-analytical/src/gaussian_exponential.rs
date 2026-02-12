use ndarray::Array1;
use statrs::function::erf::erfc;

/// Analytical convolution of exponential and standard Gaussian.
///
/// Given:
/// - f(x) = a * exp(-a * x) for x >= 0 (exponential with rate a)
/// - g(x) = exp(-0.5 * x^2) / sqrt(2 * pi) (standard Gaussian with sigma=1, mu=0)
///
/// The convolution is:
/// (f * g)(x) = 0.5 * a * exp(-x * a + 0.5 * a^2) * erfc((a - x) / sqrt(2))
pub fn gaussian_exponential_convolution(a: f64, x: f64) -> f64 {
  0.5 * a * (-x * a + 0.5 * a.powi(2)).exp() * erfc((a - x) / 2_f64.sqrt())
}

/// Evaluate gaussian-exponential convolution on grid.
pub fn gaussian_exponential_convolution_grid(a: f64, grid: &Array1<f64>) -> Array1<f64> {
  grid.mapv(|x| gaussian_exponential_convolution(a, x))
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_gaussian_exponential_convolution_positive_x() {
    let a = 0.5;
    let x = 1.0;
    let result = gaussian_exponential_convolution(a, x);
    let expected = 0.5 * a * (-x * a + 0.5 * a.powi(2)).exp() * erfc((a - x) / 2_f64.sqrt());
    assert_ulps_eq!(result, expected, max_ulps = 4);
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
