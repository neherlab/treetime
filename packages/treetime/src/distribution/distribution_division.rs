use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use crate::make_error;
use eyre::Report;
use ndarray::Array1;

const BIG_NUMBER_THRESHOLD: f64 = 1e15;
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
  let dividend_value = point.amplitude().max(TINY_NUMBER);
  let divisor_value = divisor.interp(t).unwrap_or(TINY_NUMBER).max(TINY_NUMBER);

  let dividend_log = -dividend_value.ln();
  let divisor_log = -divisor_value.ln();
  let result_log = dividend_log - divisor_log;

  if !result_log.is_finite() || result_log > BIG_NUMBER_THRESHOLD {
    return Ok(Distribution::empty());
  }

  let result_value = (-result_log).exp();
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

  let mut t_values = Vec::with_capacity(n_samples);
  let mut y_values = Vec::with_capacity(n_samples);

  let dividend_log = -(range.amplitude().max(TINY_NUMBER).ln());

  for i in 0..n_samples {
    let t = start + step * i as f64;
    let divisor_value = divisor.interp(t).unwrap_or(TINY_NUMBER).max(TINY_NUMBER);

    let divisor_log = -(divisor_value.ln());
    let result_log = dividend_log - divisor_log;

    if !result_log.is_finite() || result_log > BIG_NUMBER_THRESHOLD {
      continue;
    }

    let result_value = (-result_log).exp();
    t_values.push(t);
    y_values.push(result_value);
  }

  if t_values.is_empty() {
    return Ok(Distribution::empty());
  }

  if t_values.len() == 1 {
    return Ok(Distribution::point(t_values[0], y_values[0]));
  }

  Distribution::function(Array1::from_vec(t_values), Array1::from_vec(y_values))
}

fn divide_function_by_function(
  dividend: &DistributionFunction<f64>,
  divisor: &DistributionFunction<f64>,
) -> Result<Distribution, Report> {
  let dividend_t = dividend.t();
  let dividend_y = dividend.y();

  let n_points = dividend_t.len();
  let mut result_y = Vec::with_capacity(n_points);

  for i in 0..n_points {
    let t = dividend_t[i];
    let dividend_value = dividend_y[i].max(TINY_NUMBER);
    let divisor_value = divisor.interp(t).unwrap_or(TINY_NUMBER).max(TINY_NUMBER);

    let dividend_log = -dividend_value.ln();
    let divisor_log = -divisor_value.ln();
    let result_log = dividend_log - divisor_log;

    let result_value = if result_log.is_finite() && result_log < BIG_NUMBER_THRESHOLD {
      (-result_log).exp()
    } else {
      0.0
    };

    result_y.push(result_value);
  }

  Distribution::function(dividend_t.clone(), Array1::from_vec(result_y))
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
  fn test_divide_by_zero_fails() {
    let t = array![0.0, 1.0, 2.0];
    let y1 = array![10.0, 20.0, 30.0];
    let y2 = array![2.0, 0.0, 5.0];

    let dividend = Distribution::function(t.clone(), y1).unwrap();
    let divisor = Distribution::function(t.clone(), y2).unwrap();

    let result = distribution_division(&dividend, &divisor);
    assert!(result.is_err());
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
