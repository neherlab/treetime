use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use crate::make_error;
use eyre::Report;

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
  let range_start = range.start();
  let range_end = range.end();
  let dividend_amplitude = range.amplitude();

  let func_min = divisor.x_min();
  let func_max = divisor.x_max();

  let overlap_start = range_start.max(func_min);
  let overlap_end = range_end.min(func_max);

  if overlap_start >= overlap_end {
    return Ok(Distribution::empty());
  }

  let result_y = divisor.y().mapv(|v| dividend_amplitude / v.max(TINY_NUMBER));
  let result_fn = DistributionFunction::from_start_dx_values(func_min, divisor.dx(), result_y)?;
  Ok(Distribution::Function(result_fn))
}

fn divide_function_by_function(
  dividend: &DistributionFunction<f64>,
  divisor: &DistributionFunction<f64>,
) -> Result<Distribution, Report> {
  let div_min = dividend.x_min();
  let div_max = dividend.x_max();
  let div_dx = dividend.dx();

  let divisor_on_dividend_grid = if (divisor.x_min() - div_min).abs() < 1e-10
    && (divisor.x_max() - div_max).abs() < 1e-10
    && (divisor.dx() - div_dx).abs() < 1e-10
  {
    divisor.y().clone()
  } else {
    let resampled = divisor.resample_to_grid((div_min, div_max), div_dx)?;
    resampled.y().clone()
  };

  let safe_divisor_y = divisor_on_dividend_grid.mapv(|v| v.max(TINY_NUMBER));
  let result_y = dividend.y() / &safe_divisor_y;

  DistributionFunction::from_range_values((div_min, div_max), result_y).map(Distribution::Function)
}

#[cfg(test)]
mod tests {
  use super::*;
  use ndarray::array;
  use treetime_utils::assert_error;

  #[test]
  fn test_divide_empty_by_any() {
    let empty = Distribution::empty();
    let point = Distribution::point(1.0, 2.0);
    let actual = distribution_division(&empty, &point).unwrap();
    let expected = Distribution::empty();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_by_empty_fails() {
    let point = Distribution::point(1.0, 2.0);
    let empty = Distribution::empty();
    assert_error!(
      distribution_division(&point, &empty),
      "Cannot divide by empty distribution"
    );
  }

  #[test]
  fn test_divide_point_by_function() {
    let point = Distribution::point(2.0, 10.0);
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 5.0, 4.0, 3.0];
    let func = Distribution::function(t, y).unwrap();

    let actual = distribution_division(&point, &func).unwrap();
    let expected = Distribution::point(2.0, 2.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function() {
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y1 = array![10.0, 20.0, 30.0, 40.0, 50.0];
    let y2 = array![2.0, 4.0, 5.0, 8.0, 10.0];

    let dividend = Distribution::function(t.clone(), y1).unwrap();
    let divisor = Distribution::function(t.clone(), y2).unwrap();

    let actual = distribution_division(&dividend, &divisor).unwrap();

    let expected_y = array![5.0, 5.0, 6.0, 5.0, 5.0];
    let expected = Distribution::function(t, expected_y).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_by_zero_handled() {
    let t = array![0.0, 1.0, 2.0];
    let y1 = array![10.0, 20.0, 30.0];
    let y2 = array![2.0, 0.0, 5.0];

    let dividend = Distribution::function(t.clone(), y1).unwrap();
    let divisor = Distribution::function(t.clone(), y2).unwrap();

    let actual = distribution_division(&dividend, &divisor).unwrap();

    let expected_y = array![5.0, 20.0 / TINY_NUMBER, 6.0];
    let expected = Distribution::function(t, expected_y).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_range_by_function_full_overlap() {
    let range = Distribution::range((1.0, 3.0), 10.0);
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 5.0, 4.0, 3.0];
    let func = Distribution::function(t, y).unwrap();

    let actual = distribution_division(&range, &func).unwrap();
    let expected =
      Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![10.0, 5.0, 2.0, 2.5, 10.0 / 3.0]).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_range_by_function_partial_overlap() {
    let range = Distribution::range((1.5, 2.5), 12.0);
    let t = array![0.0, 1.0, 2.0, 3.0];
    let y = array![1.0, 2.0, 4.0, 8.0];
    let func = Distribution::function(t, y).unwrap();

    let actual = distribution_division(&range, &func).unwrap();

    let expected = Distribution::function(array![0.0, 1.0, 2.0, 3.0], array![12.0, 6.0, 3.0, 1.5]).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_range_by_function_no_overlap() {
    let range = Distribution::range((5.0, 6.0), 10.0);
    let t = array![0.0, 1.0, 2.0, 3.0];
    let y = array![1.0, 2.0, 4.0, 8.0];
    let func = Distribution::function(t, y).unwrap();

    let actual = distribution_division(&range, &func).unwrap();
    let expected = Distribution::empty();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function_same_grid() {
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y1 = array![10.0, 20.0, 30.0, 40.0, 50.0];
    let y2 = array![2.0, 4.0, 5.0, 8.0, 10.0];

    let dividend = Distribution::function(t.clone(), y1).unwrap();
    let divisor = Distribution::function(t.clone(), y2).unwrap();

    let actual = distribution_division(&dividend, &divisor).unwrap();

    let expected = Distribution::function(t, array![5.0, 5.0, 6.0, 5.0, 5.0]).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function_different_grids() {
    let t1 = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y1 = array![10.0, 20.0, 30.0, 40.0, 50.0];
    let dividend = Distribution::function(t1.clone(), y1).unwrap();

    let t2 = array![0.0, 2.0, 4.0];
    let y2 = array![2.0, 5.0, 10.0];
    let divisor = Distribution::function(t2, y2).unwrap();

    let actual = distribution_division(&dividend, &divisor).unwrap();

    let expected = Distribution::function(t1, array![5.0, 20.0 / 3.5, 6.0, 40.0 / 7.5, 5.0]).unwrap();
    assert_eq!(expected, actual);
  }
}
