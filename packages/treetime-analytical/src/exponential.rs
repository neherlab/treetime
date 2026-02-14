use ndarray::Array1;

/// Exponential PDF: f(x) = rate * exp(-rate * x) for x >= 0, 0 otherwise.
///
/// # Preconditions
///
/// `rate` must be positive. Non-positive rates do not define a valid exponential distribution.
pub fn exponential_pdf(rate: f64, x: f64) -> f64 {
  debug_assert!(rate > 0.0, "exponential_pdf: rate must be positive");
  if x < 0.0 { 0.0 } else { rate * (-rate * x).exp() }
}

/// Evaluate exponential PDF on grid.
pub fn exponential_pdf_grid(rate: f64, grid: &Array1<f64>) -> Array1<f64> {
  grid.mapv(|x| exponential_pdf(rate, x))
}

/// Analytical convolution of two exponential distributions.
///
/// Given f(x) = a * exp(-a * x) and g(x) = b * exp(-b * x) for x >= 0,
/// the convolution (f * g)(x) is:
/// - For a != b: (a * b) / (a - b) * (1 - exp(-(a - b) * x)) * exp(-b * x)
/// - For a = b (limit form): a * b * x * exp(-a * x)
///
/// # Preconditions
///
/// Both `a` and `b` must be positive. Non-positive rates do not define valid exponential distributions.
pub fn exponential_convolution(a: f64, b: f64, x: f64) -> f64 {
  debug_assert!(a > 0.0 && b > 0.0, "exponential_convolution: rates must be positive");
  if x < 0.0 {
    0.0
  } else if (a - b).abs() < 1e-15 {
    a * b * x * (-a * x).exp()
  } else {
    // Use exp_m1 to avoid catastrophic cancellation when a ≈ b.
    // Since 1 - exp(y) = -exp_m1(y), we have:
    // (1 - exp(-(a-b)*x)) = -exp_m1(-(a-b)*x)
    // Formula: (a*b)/(a-b) * (1 - exp(-(a-b)*x)) * exp(-b*x)
    //        = (a*b)/(b-a) * exp_m1(-(a-b)*x) * exp(-b*x)
    (a * b) / (b - a) * (-(a - b) * x).exp_m1() * (-b * x).exp()
  }
}

/// Evaluate exponential convolution on grid.
pub fn exponential_convolution_grid(a: f64, b: f64, grid: &Array1<f64>) -> Array1<f64> {
  grid.mapv(|x| exponential_convolution(a, b, x))
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_exponential_pdf_positive() {
    let rate = 2.0;
    let result = exponential_pdf(rate, 1.0);
    let expected = 2.0 * (-2.0_f64).exp();
    assert_ulps_eq!(result, expected, max_ulps = 4);
  }

  #[test]
  fn test_exponential_pdf_negative() {
    let result = exponential_pdf(1.0, -1.0);
    assert_ulps_eq!(result, 0.0, max_ulps = 4);
  }

  #[test]
  fn test_exponential_pdf_zero() {
    let rate = 2.0;
    let result = exponential_pdf(rate, 0.0);
    assert_ulps_eq!(result, rate, max_ulps = 4);
  }

  #[test]
  fn test_exponential_convolution_distinct_rates() {
    let a = 1.0;
    let b = 2.0;
    let x = 1.0;
    let result = exponential_convolution(a, b, x);
    let expected = (a * b) / (a - b) * (1.0 - (-(a - b) * x).exp()) * (-b * x).exp();
    assert_ulps_eq!(result, expected, max_ulps = 4);
  }

  #[test]
  fn test_exponential_convolution_equal_rates() {
    let a = 1.5;
    let b = 1.5;
    let x = 2.0;
    let result = exponential_convolution(a, b, x);
    let expected = a * b * x * (-a * x).exp();
    assert_ulps_eq!(result, expected, max_ulps = 4);
  }

  #[test]
  fn test_exponential_convolution_negative_x() {
    let result = exponential_convolution(1.0, 2.0, -1.0);
    assert_ulps_eq!(result, 0.0, max_ulps = 4);
  }

  #[test]
  fn test_exponential_convolution_at_zero() {
    let result = exponential_convolution(1.0, 2.0, 0.0);
    assert_ulps_eq!(result, 0.0, max_ulps = 4);
  }

  #[test]
  fn test_exponential_pdf_grid() {
    let grid = array![-1.0, 0.0, 1.0, 2.0];
    let rate = 1.0;
    let result = exponential_pdf_grid(rate, &grid);

    assert_ulps_eq!(result[0], 0.0, max_ulps = 4);
    assert_ulps_eq!(result[1], 1.0, max_ulps = 4);
    assert_ulps_eq!(result[2], (-1.0_f64).exp(), max_ulps = 4);
    assert_ulps_eq!(result[3], (-2.0_f64).exp(), max_ulps = 4);
  }

  #[test]
  fn test_exponential_convolution_grid() {
    let grid = array![-1.0, 0.0, 1.0];
    let result = exponential_convolution_grid(1.0, 2.0, &grid);

    assert_ulps_eq!(result[0], 0.0, max_ulps = 4);
    assert_ulps_eq!(result[1], 0.0, max_ulps = 4);
    assert!(result[2] > 0.0);
  }

  #[test]
  fn test_exponential_convolution_near_equal_rates() {
    // Test numerical stability when rates are nearly equal.
    // The exp_m1 formulation should avoid catastrophic cancellation.
    let a = 1.0;
    let b = 1.0 + 1e-10;
    let x = 2.0;
    let result = exponential_convolution(a, b, x);

    // Result should be close to the equal-rates limit: a * b * x * exp(-a * x)
    let limit = a * b * x * (-a * x).exp();
    assert!((result - limit).abs() < 1e-6, "result={result}, limit={limit}");
    assert!(result.is_finite(), "result should be finite");
  }
}
