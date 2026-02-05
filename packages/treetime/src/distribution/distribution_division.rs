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
  let values_array: ndarray::Array1<f64> = filtered.iter().map(|(_, v)| *v).collect();

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
