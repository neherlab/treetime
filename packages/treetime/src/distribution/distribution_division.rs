use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use crate::distribution::y_axis_policy::YAxisPolicy;
use crate::make_error;
use eyre::Report;

pub fn distribution_division<Y: YAxisPolicy>(
  dividend: &Distribution<Y>,
  divisor: &Distribution<Y>,
) -> Result<Distribution<Y>, Report> {
  match (dividend, divisor) {
    (Distribution::Empty, _) => {
      Ok(Distribution::Empty) //
    },
    (_, Distribution::Empty) => {
      make_error!("Cannot divide by empty distribution") //
    },
    (Distribution::Point(_), Distribution::Point(_)) => {
      make_error!("Cannot divide Point by Point: operation not well-defined") //
    },
    (Distribution::Range(_), Distribution::Point(_)) => {
      make_error!("Cannot divide Range by Point: operation not well-defined") //
    },
    (Distribution::Function(_), Distribution::Point(_)) => {
      make_error!("Cannot divide Function by Point: operation not well-defined") //
    },
    (Distribution::Point(_), Distribution::Range(_)) => {
      make_error!("Cannot divide Point by Range: operation not well-defined") //
    },
    (Distribution::Range(_), Distribution::Range(_)) => {
      make_error!("Cannot divide Range by Range: operation not well-defined") //
    },
    (Distribution::Function(_), Distribution::Range(_)) => {
      make_error!("Cannot divide Function by Range: operation not well-defined") //
    },
    (Distribution::Point(a), Distribution::Function(b)) => {
      divide_point_by_function::<Y>(a, b) //
    },
    (Distribution::Range(a), Distribution::Function(b)) => {
      divide_range_by_function::<Y>(a, b) //
    },
    (Distribution::Function(a), Distribution::Function(b)) => {
      divide_function_by_function::<Y>(a, b) //
    },
    (Distribution::Formula(_), _) | (_, Distribution::Formula(_)) => {
      panic!("Division not implemented for Formula distributions")
    }, //
  }
}

fn divide_point_by_function<Y: YAxisPolicy>(
  point: &DistributionPoint<f64, Y>,
  divisor: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let t = point.t();
  let dividend_value = point.amplitude();
  let divisor_value = divisor.interp(t).unwrap_or_else(|_| Y::multiplicative_identity());
  let safe_divisor_value = Y::safe_divisor(divisor_value);
  let result_value = Y::divide(dividend_value, safe_divisor_value);

  if !Y::is_defined(result_value) {
    return Ok(Distribution::empty());
  }

  Ok(Distribution::point(t, result_value))
}

fn divide_range_by_function<Y: YAxisPolicy>(
  range: &DistributionRange<f64, Y>,
  divisor: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  const EPS: f64 = 1e-9;

  let dividend_amplitude = range.amplitude();
  let func_t = divisor.t();
  let func_y = divisor.y();

  // Filter to range support (consistent with multiply_range_function)
  let filtered: Vec<(f64, f64)> = func_t
    .iter()
    .zip(func_y.iter())
    .filter(|&(&t, _)| t >= range.start() - EPS && t <= range.end() + EPS)
    .map(|(&t, &y)| (t, Y::divide(dividend_amplitude, Y::safe_divisor(y))))
    .collect();

  if filtered.is_empty() {
    return Ok(Distribution::empty());
  }

  // Single point -> return Point distribution
  if filtered.len() == 1 {
    let (t, v) = filtered[0];
    return Ok(Distribution::point(t, v));
  }

  let overlap_min = filtered.first().unwrap().0;
  let overlap_max = filtered.last().unwrap().0;
  let values: Vec<f64> = filtered.iter().map(|(_, v)| *v).collect();
  let values_array = ndarray::Array1::from_vec(values);

  let result_fn = DistributionFunction::from_range_values((overlap_min, overlap_max), values_array)?;
  Ok(Distribution::Function(result_fn))
}

fn divide_function_by_function<Y: YAxisPolicy>(
  dividend: &DistributionFunction<f64, Y>,
  divisor: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let div_min = dividend.x_min();
  let div_max = dividend.x_max();
  let div_dx = dividend.dx();

  let divisor_on_dividend_grid = if (divisor.x_min() - div_min).abs() < 1e-10
    && (divisor.x_max() - div_max).abs() < 1e-10
    && (divisor.dx() - div_dx).abs() < 1e-10
  {
    divisor.y().clone()
  } else {
    let resampled = divisor.resample_range_dx((div_min, div_max), div_dx)?;
    resampled.y().clone()
  };

  let result_y = dividend
    .y()
    .iter()
    .zip(divisor_on_dividend_grid.iter())
    .map(|(&d, &s)| Y::divide(d, Y::safe_divisor(s)))
    .collect();

  DistributionFunction::from_range_values((div_min, div_max), result_y).map(Distribution::Function)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::distribution::distribution::DistributionPlain as Distribution;
  use ndarray::array;
  use treetime_utils::assert_error;

  const TINY_NUMBER: f64 = 1e-10;

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
    // Range [1.0, 3.0] overlaps with function at points [1.0, 2.0, 3.0]
    let range = Distribution::range((1.0, 3.0), 10.0);
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 5.0, 4.0, 3.0];
    let func = Distribution::function(t, y).unwrap();

    let actual = distribution_division(&range, &func).unwrap();
    // Result should only cover the overlap region [1.0, 3.0]
    let expected = Distribution::function(array![1.0, 2.0, 3.0], array![5.0, 2.0, 2.5]).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_range_by_function_partial_overlap() {
    // Range [1.5, 2.5] overlaps with function at point [2.0]
    let range = Distribution::range((1.5, 2.5), 12.0);
    let t = array![0.0, 1.0, 2.0, 3.0];
    let y = array![1.0, 2.0, 4.0, 8.0];
    let func = Distribution::function(t, y).unwrap();

    let actual = distribution_division(&range, &func).unwrap();
    // Result should only cover the overlap region (point at 2.0)
    let expected = Distribution::point(2.0, 3.0);
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
