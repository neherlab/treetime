use crate::Distribution;
use crate::distribution_core::function::DistributionFunction;
use crate::distribution_core::point::DistributionPoint;
use crate::distribution_core::range::DistributionRange;
use crate::distribution_ops::multiply::{
  HardDomain, distribution_hard_domain, function_hard_domain, guarded_empty_result,
  multiplication_support_intersection, point_hard_domain, range_hard_domain,
};
use crate::distribution_ops::time_bounds::{SupportIntersection, distribution_support_n_points};
use crate::policy::YAxisPolicy;
use eyre::Report;
use ndarray::{Array1, Zip};
use std::cmp::Ordering;
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, HardApproachLaw, HardCell, MaxGridPoints, Side};
use treetime_utils::make_error;

pub fn distribution_division<Y: YAxisPolicy>(
  dividend: &Distribution<Y>,
  divisor: &Distribution<Y>,
  max_points: MaxGridPoints,
) -> Result<Distribution<Y>, Report> {
  match (dividend, divisor) {
    (Distribution::Formula(_), _) | (_, Distribution::Formula(_)) => {
      make_error!("Cannot divide {dividend} by {divisor}: operation not implemented")
    },
    (Distribution::Empty, _) => {
      let dividend_domain = distribution_hard_domain(dividend);
      let divisor_domain = distribution_hard_domain(divisor);
      guarded_empty_result("division", dividend_domain, divisor_domain)
    },
    (_, Distribution::Empty) => make_error!("Cannot divide by empty distribution"),
    (
      Distribution::Point(_) | Distribution::Range(_) | Distribution::Function(_),
      Distribution::Point(_) | Distribution::Range(_),
    ) => make_error!("Cannot divide {dividend} by {divisor}: operation not well-defined"),
    (Distribution::Point(a), Distribution::Function(b)) => divide_point_by_function::<Y>(a, b),
    (Distribution::Range(a), Distribution::Function(b)) => divide_range_by_function::<Y>(a, b, max_points),
    (Distribution::Function(a), Distribution::Function(b)) => divide_function_by_function::<Y>(a, b, max_points),
  }
}

fn divide_point_by_function<Y: YAxisPolicy>(
  point: &DistributionPoint<f64, Y>,
  divisor: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let t = point.t();
  let tail = if t < divisor.x_min() {
    Some(divisor.left_extrap())
  } else if t > divisor.x_max() {
    Some(divisor.right_extrap())
  } else {
    None
  };
  if let Some(tail) = tail {
    if !tail.is_soft() {
      return guarded_empty_result(
        "division",
        Some(point_hard_domain(point)),
        Some(function_hard_domain(divisor)),
      );
    }
  }

  let divisor_value = divisor.interp(t)?;
  let result_value = Y::divide(point.amplitude(), Y::safe_divisor(divisor_value));
  if !Y::is_defined(result_value) {
    return guarded_empty_result(
      "division",
      Some(point_hard_domain(point)),
      Some(function_hard_domain(divisor)),
    );
  }

  Ok(Distribution::point(t, result_value))
}

fn divide_range_by_function<Y: YAxisPolicy>(
  range: &DistributionRange<f64, Y>,
  divisor: &DistributionFunction<f64, Y>,
  max_points: MaxGridPoints,
) -> Result<Distribution<Y>, Report> {
  let range_domain = range_hard_domain(range);
  let divisor_domain = function_hard_domain(divisor);
  match multiplication_support_intersection(&[range_domain, divisor_domain]) {
    SupportIntersection::Disjoint => {
      let dividend_domain = Some(range_domain);
      let divisor_domain = Some(divisor_domain);
      guarded_empty_result("division", dividend_domain, divisor_domain)
    },
    SupportIntersection::Point(t) => Ok(Distribution::point(
      t,
      Y::divide(range.amplitude(), Y::safe_divisor(divisor.interp(t)?)),
    )),
    SupportIntersection::Interval(bounds) => {
      let n_points = distribution_support_n_points(bounds, divisor.dx(), max_points)?;
      let grid = Array1::linspace(bounds.0, bounds.1, n_points);
      let values = divisor
        .interp_many(&grid)?
        .mapv(|value| Y::divide(range.amplitude(), Y::safe_divisor(value)));
      let function = DistributionFunction::from_range_values(bounds, values)?;
      let function = apply_division_tail(function, Side::Left, range_domain, divisor_domain)?;
      let function = apply_division_tail(function, Side::Right, range_domain, divisor_domain)?;
      Ok(Distribution::Function(function))
    },
  }
}

fn divide_function_by_function<Y: YAxisPolicy>(
  dividend: &DistributionFunction<f64, Y>,
  divisor: &DistributionFunction<f64, Y>,
  max_points: MaxGridPoints,
) -> Result<Distribution<Y>, Report> {
  let dividend_domain = function_hard_domain(dividend);
  let divisor_domain = function_hard_domain(divisor);
  match multiplication_support_intersection(&[dividend_domain, divisor_domain]) {
    SupportIntersection::Disjoint => {
      let dividend_domain = Some(dividend_domain);
      let divisor_domain = Some(divisor_domain);
      guarded_empty_result("division", dividend_domain, divisor_domain)
    },
    SupportIntersection::Point(t) => Ok(Distribution::point(
      t,
      Y::divide(dividend.interp(t)?, Y::safe_divisor(divisor.interp(t)?)),
    )),
    SupportIntersection::Interval(bounds) => {
      let n_points = distribution_support_n_points(bounds, dividend.dx().min(divisor.dx()), max_points)?;
      let grid = Array1::linspace(bounds.0, bounds.1, n_points);
      let dividend_values = dividend.interp_many_stored(&grid)?;
      let divisor_values = divisor.interp_many_stored(&grid)?;
      let mut values = Zip::from(&dividend_values)
        .and(&divisor_values)
        .map_collect(|&dividend, &divisor| Y::divide(dividend, Y::safe_divisor(divisor)));
      let last = values.len() - 1;
      let left_end = shared_hard_end::<Y>(dividend, divisor, Side::Left, bounds.0)?;
      let right_end = shared_hard_end::<Y>(dividend, divisor, Side::Right, bounds.1)?;
      if let Some(end) = left_end {
        values[0] = end.value::<Y>(values[1]);
      }
      if let Some(end) = right_end {
        values[last] = end.value::<Y>(values[last - 1]);
      }
      let function = DistributionFunction::from_range_values(bounds, values)?
        .with_amplitude(Y::divide(dividend.amplitude(), divisor.amplitude()));
      let function = apply_division_tail(function, Side::Left, dividend_domain, divisor_domain)?;
      let function = apply_division_tail(function, Side::Right, dividend_domain, divisor_domain)?;
      let function = match left_end.and_then(SharedHardEnd::tail) {
        Some(tail) => function.with_left_extrap(tail)?,
        None => function,
      };
      let function = match right_end.and_then(SharedHardEnd::tail) {
        Some(tail) => function.with_right_extrap(tail)?,
        None => function,
      };
      Ok(Distribution::Function(function))
    },
  }
}

fn shared_hard_end<Y: YAxisPolicy>(
  dividend: &DistributionFunction<f64, Y>,
  divisor: &DistributionFunction<f64, Y>,
  side: Side,
  t: f64,
) -> Result<Option<SharedHardEnd>, Report> {
  if divisor.interp_stored(t)?.is_finite() {
    return Ok(None);
  }
  let (Some(dividend_law), Some(divisor_law)) = (end_approach(dividend, side, t)?, end_approach(divisor, side, t)?)
  else {
    return Ok(None);
  };
  let exponent = dividend_law.exponent - divisor_law.exponent;
  let end = match exponent.total_cmp(&0.0) {
    Ordering::Greater => SharedHardEnd::ZeroDensity(HardApproachLaw { t_hard: t, b: exponent }),
    Ordering::Equal => SharedHardEnd::Finite(Y::divide(dividend_law.offset, divisor_law.offset)),
    Ordering::Less => SharedHardEnd::Flat,
  };
  Ok(Some(end))
}

fn end_approach<Y: YAxisPolicy>(
  f: &DistributionFunction<f64, Y>,
  side: Side,
  t: f64,
) -> Result<Option<EndApproach>, Report> {
  let value = f.interp_stored(t)?;
  if value.is_finite() {
    return Ok(Some(EndApproach {
      exponent: 0.0,
      offset: value,
    }));
  }
  let Some(last_cell) = f.len().checked_sub(2) else {
    return Ok(None);
  };
  let (edge, cell) = match side {
    Side::Left => (f.x_min(), 0),
    Side::Right => (f.x_max(), last_cell),
  };
  if !f.grid_fn().coincides(edge, t) {
    return Ok(None);
  }
  let Some(HardCell { law, anchor }) = f.grid_fn().hard_cell(cell, f.left_extrap(), f.right_extrap())? else {
    return Ok(None);
  };
  Ok(Some(EndApproach {
    exponent: law.b,
    offset: anchor.y + law.b * (anchor.t - law.t_hard).abs().ln(),
  }))
}

#[derive(Clone, Copy, Debug)]
struct EndApproach {
  exponent: f64,
  offset: f64,
}

#[derive(Clone, Copy, Debug)]
enum SharedHardEnd {
  ZeroDensity(HardApproachLaw),
  Finite(f64),
  Flat,
}

impl SharedHardEnd {
  fn value<Y: YAxisPolicy>(self, inner: f64) -> f64 {
    match self {
      SharedHardEnd::ZeroDensity(_) => Y::probability_zero(),
      SharedHardEnd::Finite(value) => value,
      SharedHardEnd::Flat => inner,
    }
  }

  fn tail(self) -> Option<BoundaryBehavior> {
    match self {
      SharedHardEnd::ZeroDensity(law) => Some(BoundaryBehavior::HardApproach(law)),
      SharedHardEnd::Finite(_) | SharedHardEnd::Flat => None,
    }
  }
}

fn apply_division_tail<Y: YAxisPolicy>(
  function: DistributionFunction<f64, Y>,
  side: Side,
  dividend: HardDomain,
  divisor: HardDomain,
) -> Result<DistributionFunction<f64, Y>, Report> {
  match division_side_tail(side, dividend, divisor) {
    DivisionSideTail::Refit => function.fit_soft_tail(side, DEFAULT_TAIL_FIT_POINTS),
    DivisionSideTail::Fixed(tail) => match side {
      Side::Left => function.with_left_extrap(tail),
      Side::Right => function.with_right_extrap(tail),
    },
  }
}

fn division_side_tail(side: Side, dividend: HardDomain, divisor: HardDomain) -> DivisionSideTail {
  let dividend_tail = side_tail(side, dividend);
  let divisor_tail = side_tail(side, divisor);
  match (dividend_tail.is_soft(), divisor_tail.is_soft()) {
    (true, true) => DivisionSideTail::Refit,
    (true, false) => DivisionSideTail::Fixed(BoundaryBehavior::Error),
    (false, true) => DivisionSideTail::Fixed(dividend_tail),
    (false, false) => {
      let dividend_inner = match side {
        Side::Left => side_bound(side, dividend) >= side_bound(side, divisor),
        Side::Right => side_bound(side, dividend) <= side_bound(side, divisor),
      };
      if dividend_inner {
        DivisionSideTail::Fixed(dividend_tail)
      } else {
        DivisionSideTail::Fixed(BoundaryBehavior::Error)
      }
    },
  }
}

fn side_tail(side: Side, domain: HardDomain) -> BoundaryBehavior {
  match side {
    Side::Left => domain.1.0,
    Side::Right => domain.1.1,
  }
}

fn side_bound(side: Side, domain: HardDomain) -> f64 {
  match side {
    Side::Left => domain.0.0,
    Side::Right => domain.0.1,
  }
}

enum DivisionSideTail {
  Refit,
  Fixed(BoundaryBehavior),
}
