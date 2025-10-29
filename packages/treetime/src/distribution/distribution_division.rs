use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use crate::make_error;
use eyre::Report;
use ndarray::Array1;

const TINY_NUMBER: f64 = 1e-10;

pub fn distribution_division(dividend: &Distribution, divisor: &Distribution) -> Result<Distribution, Report> {
  match (dividend, divisor) {
    (Distribution::Empty, _) => Ok(Distribution::Empty),
    (_, Distribution::Empty) => make_error!("Cannot divide by empty distribution"),
    (Distribution::Point(_), Distribution::Point(_)) => {
      make_error!("Cannot divide Point by Point: operation not well-defined")
    },
    (Distribution::Range(_), Distribution::Point(_)) => {
      make_error!("Cannot divide Range by Point: operation not well-defined")
    },
    (Distribution::Function(_), Distribution::Point(_)) => {
      make_error!("Cannot divide Function by Point: operation not well-defined")
    },
    (Distribution::Point(_), Distribution::Range(_)) => {
      make_error!("Cannot divide Point by Range: operation not well-defined")
    },
    (Distribution::Range(_), Distribution::Range(_)) => {
      make_error!("Cannot divide Range by Range: operation not well-defined")
    },
    (Distribution::Function(_), Distribution::Range(_)) => {
      make_error!("Cannot divide Function by Range: operation not well-defined")
    },
    (Distribution::Point(a), Distribution::Function(b)) => divide_point_by_function(a, b),
    (Distribution::Range(a), Distribution::Function(b)) => divide_range_by_function(a, b),
    (Distribution::Function(a), Distribution::Function(b)) => divide_function_by_function(a, b),
  }
}

fn divide_point_by_function(
  point: &DistributionPoint<f64>,
  divisor: &DistributionFunction<f64>,
) -> Result<Distribution, Report> {
  let t = point.t();
  let dividend_value = point.amplitude();
  let divisor_value = divisor.interp(t).unwrap_or(TINY_NUMBER).max(TINY_NUMBER);

  let result_value = dividend_value / divisor_value;

  if !result_value.is_finite() {
    return Ok(Distribution::empty());
  }

  Ok(Distribution::point(t, result_value))
}

fn divide_range_by_function(
  range: &DistributionRange<f64>,
  divisor: &DistributionFunction<f64>,
) -> Result<Distribution, Report> {
  let n_samples = 100;
  let start = range.start();
  let end = range.end();
  let step = (end - start) / (n_samples - 1) as f64;

  let dividend_amplitude = range.amplitude();

  // Use ndarray to generate uniform t values
  let t_values = Array1::range(0.0, n_samples as f64, 1.0)
    .mapv(|i| start + step * i);

  // Compute divisor values using interp_many for efficiency
  let divisor_values = divisor.interp_many(&t_values)?;

  // Apply TINY_NUMBER safety and perform division
  let safe_divisor_values = divisor_values.mapv(|v| v.max(TINY_NUMBER));
  let result_y = Array1::from_elem(n_samples, dividend_amplitude) / &safe_divisor_values;

  Distribution::function(t_values, result_y)
}

fn divide_function_by_function(
  dividend: &DistributionFunction<f64>,
  divisor: &DistributionFunction<f64>,
) -> Result<Distribution, Report> {
  let dividend_t = dividend.t();
  let dividend_y = dividend.y();

  // Use ndarray for vectorized operations
  let divisor_values = divisor.interp_many(dividend_t)?;

  // Apply TINY_NUMBER safety for zero or near-zero values
  let safe_divisor_values = divisor_values.mapv(|v| v.max(TINY_NUMBER));

  // Perform vectorized division
  let result_y = dividend_y / &safe_divisor_values;

  Distribution::function(dividend_t.clone(), result_y)
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_divide_empty_by_any() {
    let empty = Distribution::empty();
    let point = Distribution::point(1.0, 2.0);
    let result = distribution_division(&empty, &point).unwrap();
    assert_eq!(result, Distribution::empty());
  }

  #[test]
  fn test_divide_by_empty_fails() {
    let point = Distribution::point(1.0, 2.0);
    let empty = Distribution::empty();
    let result = distribution_division(&point, &empty);
    assert!(result.is_err());
  }

  #[test]
  fn test_divide_point_by_function() {
    let point = Distribution::point(2.0, 10.0);
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 5.0, 4.0, 3.0];
    let func = Distribution::function(t, y).unwrap();

    let result = distribution_division(&point, &func).unwrap();

    match result {
      Distribution::Point(p) => {
        assert_ulps_eq!(p.t(), 2.0);
        assert_ulps_eq!(p.amplitude(), 2.0);
      },
      _ => panic!("Expected Point distribution"),
    }
  }

  #[test]
  fn test_divide_function_by_function() {
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y1 = array![10.0, 20.0, 30.0, 40.0, 50.0];
    let y2 = array![2.0, 4.0, 5.0, 8.0, 10.0];

    let dividend = Distribution::function(t.clone(), y1).unwrap();
    let divisor = Distribution::function(t.clone(), y2).unwrap();

    let result = distribution_division(&dividend, &divisor).unwrap();

    match result {
      Distribution::Function(f) => {
        assert_ulps_eq!(f.y()[0], 5.0);
        assert_ulps_eq!(f.y()[1], 5.0);
        assert_ulps_eq!(f.y()[2], 6.0);
        assert_ulps_eq!(f.y()[3], 5.0);
        assert_ulps_eq!(f.y()[4], 5.0);
      },
      _ => panic!("Expected Function distribution"),
    }
  }

  #[test]
  fn test_divide_by_zero_handled() {
    let t = array![0.0, 1.0, 2.0];
    let y1 = array![10.0, 20.0, 30.0];
    let y2 = array![2.0, 0.0, 5.0];

    let dividend = Distribution::function(t.clone(), y1).unwrap();
    let divisor = Distribution::function(t.clone(), y2).unwrap();

    let result = distribution_division(&dividend, &divisor).unwrap();
    // Should succeed with TINY_NUMBER handling for zero divisor, maintaining uniform grid
    match result {
      Distribution::Function(f) => {
        assert_eq!(f.t().len(), 3);
        assert_ulps_eq!(f.y()[0], 5.0); // 10.0 / 2.0
        assert!(f.y()[1] > 1e9); // 20.0 / TINY_NUMBER (very large)
        assert_ulps_eq!(f.y()[2], 6.0); // 30.0 / 5.0
      },
      _ => panic!("Expected Function distribution"),
    }
  }

  #[test]
  fn test_divide_range_by_function() {
    let range = Distribution::range((1.0, 3.0), 10.0);
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 5.0, 4.0, 3.0];
    let func = Distribution::function(t, y).unwrap();

    let result = distribution_division(&range, &func).unwrap();

    match result {
      Distribution::Function(f) => {
        assert!(f.t().len() > 0);
        assert!(f.y().len() > 0);
        assert_eq!(f.t().len(), f.y().len());
      },
      _ => panic!("Expected Function distribution"),
    }
  }
}
