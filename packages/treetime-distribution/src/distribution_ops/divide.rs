use crate::Distribution;
use crate::distribution_core::function::DistributionFunction;
use crate::distribution_core::point::DistributionPoint;
use crate::distribution_core::range::DistributionRange;
use crate::distribution_ops::multiply::{
  HardDomain, distribution_hard_domain, guarded_empty_result, point_hard_domain, range_hard_domain,
};
use crate::distribution_ops::time_bounds::{
  SupportIntersection, distribution_support_intersection, distribution_support_n_points,
};
use crate::policy::YAxisPolicy;
use eyre::Report;
use ndarray::{Array1, Zip};
use treetime_grid::BoundaryBehavior;
use treetime_utils::make_error;

pub fn distribution_division<Y: YAxisPolicy>(
  dividend: &Distribution<Y>,
  divisor: &Distribution<Y>,
) -> Result<Distribution<Y>, Report> {
  match (dividend, divisor) {
    (Distribution::Formula(_), Distribution::Empty) => {
      make_error!("Cannot divide Formula by Empty: operation not implemented")
    },
    (Distribution::Formula(_), Distribution::Point(_)) => {
      make_error!("Cannot divide Formula by Point: operation not implemented")
    },
    (Distribution::Formula(_), Distribution::Range(_)) => {
      make_error!("Cannot divide Formula by Range: operation not implemented")
    },
    (Distribution::Formula(_), Distribution::Function(_)) => {
      make_error!("Cannot divide Formula by Function: operation not implemented")
    },
    (Distribution::Formula(_), Distribution::Formula(_)) => {
      make_error!("Cannot divide Formula by Formula: operation not implemented")
    },
    (Distribution::Empty, Distribution::Formula(_)) => {
      make_error!("Cannot divide Empty by Formula: operation not implemented")
    },
    (Distribution::Point(_), Distribution::Formula(_)) => {
      make_error!("Cannot divide Point by Formula: operation not implemented")
    },
    (Distribution::Range(_), Distribution::Formula(_)) => {
      make_error!("Cannot divide Range by Formula: operation not implemented")
    },
    (Distribution::Function(_), Distribution::Formula(_)) => {
      make_error!("Cannot divide Function by Formula: operation not implemented")
    },
    (Distribution::Empty, _) => guarded_empty_result(
      "division",
      distribution_hard_domain(dividend),
      distribution_hard_domain(divisor),
    ),
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
  }
}

fn divide_point_by_function<Y: YAxisPolicy>(
  point: &DistributionPoint<f64, Y>,
  divisor: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let t = point.t();
  let dividend_value = point.amplitude();
  let divisor_value = divisor.interp(t)?;
  let safe_divisor_value = Y::safe_divisor(divisor_value);
  let result_value = Y::divide(dividend_value, safe_divisor_value);

  if !Y::is_defined(result_value) {
    return guarded_empty_result(
      "division",
      Some(point_hard_domain(point)),
      Some(division_divisor_domain(divisor)),
    );
  }

  Ok(Distribution::point(t, result_value))
}

fn divide_range_by_function<Y: YAxisPolicy>(
  range: &DistributionRange<f64, Y>,
  divisor: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  match division_support_intersection((range.start(), range.end()), divisor) {
    SupportIntersection::Disjoint => guarded_empty_result(
      "division",
      Some(range_hard_domain(range)),
      Some(division_divisor_domain(divisor)),
    ),
    SupportIntersection::Point(t) => Ok(Distribution::point(
      t,
      Y::divide(range.amplitude(), Y::safe_divisor(divisor.interp(t)?)),
    )),
    SupportIntersection::Interval(bounds) => {
      let n_points = distribution_support_n_points(bounds, divisor.dx())?;
      let grid = Array1::linspace(bounds.0, bounds.1, n_points);
      let values = divisor
        .interp_many(&grid)?
        .mapv(|value| Y::divide(range.amplitude(), Y::safe_divisor(value)));
      let function = DistributionFunction::from_range_values(bounds, values)?;
      Ok(Distribution::Function(function))
    },
  }
}

fn divide_function_by_function<Y: YAxisPolicy>(
  dividend: &DistributionFunction<f64, Y>,
  divisor: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  match division_support_intersection((dividend.x_min(), dividend.x_max()), divisor) {
    SupportIntersection::Disjoint => guarded_empty_result(
      "division",
      Some((
        (dividend.x_min(), dividend.x_max()),
        (BoundaryBehavior::Error, BoundaryBehavior::Error),
      )),
      Some(division_divisor_domain(divisor)),
    ),
    SupportIntersection::Point(t) => Ok(Distribution::point(
      t,
      Y::divide(dividend.interp(t)?, Y::safe_divisor(divisor.interp(t)?)),
    )),
    SupportIntersection::Interval(bounds) => {
      let n_points = distribution_support_n_points(bounds, dividend.dx().min(divisor.dx()))?;
      let grid = Array1::linspace(bounds.0, bounds.1, n_points);
      let dividend_values = dividend.interp_many(&grid)?;
      let divisor_values = divisor.interp_many(&grid)?;
      let values = Zip::from(&dividend_values)
        .and(&divisor_values)
        .map_collect(|&dividend, &divisor| Y::divide(dividend, Y::safe_divisor(divisor)));
      // Tail policy survives (kb/decisions/distribution-tails-and-arithmetic.md). Division is
      // asymmetric: the dividend's bounds are used as-is, and every divisor tail bounds the quotient
      // at the divisor's grid edge, so the divisor is never sampled past its own support. The
      // quotient inherits the dividend's tail: beyond a soft dividend edge the quotient keeps
      // decaying, beyond a hard edge it is zero. A divisor `Error` edge strictly inside the dividend
      // truncates the result grid there, and the quotient is undefined past it, so that side
      // becomes `Error`.
      let left_truncated = divisor.left_extrap() == BoundaryBehavior::Error && divisor.x_min() > dividend.x_min();
      let right_truncated = divisor.right_extrap() == BoundaryBehavior::Error && divisor.x_max() < dividend.x_max();
      let function = DistributionFunction::from_range_values(bounds, values)?
        .with_left_extrap(division_result_tail(dividend.left_extrap(), left_truncated))?
        .with_right_extrap(division_result_tail(dividend.right_extrap(), right_truncated))?;
      Ok(Distribution::Function(function))
    },
  }
}

/// Per-side tail of a Function/Function quotient: the dividend's tail, unless a divisor `Error` bound
/// truncated the result grid inside the dividend on that side (then the quotient is undefined beyond
/// the truncating edge, so `Error`).
fn division_result_tail(dividend_tail: BoundaryBehavior, divisor_truncates: bool) -> BoundaryBehavior {
  if divisor_truncates {
    BoundaryBehavior::Error
  } else {
    dividend_tail
  }
}

fn division_support_intersection<Y: YAxisPolicy>(
  dividend_bounds: (f64, f64),
  divisor: &DistributionFunction<f64, Y>,
) -> SupportIntersection {
  // Every divisor tail bounds the quotient at the divisor's real, grid-backed edge, so the divisor
  // is never sampled past its own support. Dividing by an extrapolated divisor tail beyond the grid
  // would produce a spurious spike: a `Hard`/`HardApproach` edge is zero probability (`+inf` under
  // `NegLog`, giving `dividend - (+inf) = -inf`); a `Linear` edge decays, so its reciprocal
  // inflates the quotient; an `Error` edge is undefined. The dividend's own tails carry the cavity
  // beyond the divisor's support. See `kb/decisions/distribution-tails-and-arithmetic.md`
  // (Division).
  let divisor_bounds = (divisor.x_min(), divisor.x_max());
  distribution_support_intersection(dividend_bounds, divisor_bounds)
}

/// The divisor's hard domain as division treats it, so [`guarded_empty_result`] reproduces exactly
/// the disjointness decided by [`division_support_intersection`]. Every divisor tail bounds the
/// quotient at the divisor's grid edge (never sampled past its support), so both sides are hard at
/// the grid bounds.
fn division_divisor_domain<Y: YAxisPolicy>(divisor: &DistributionFunction<f64, Y>) -> HardDomain {
  (
    (divisor.x_min(), divisor.x_max()),
    (BoundaryBehavior::Error, BoundaryBehavior::Error),
  )
}
