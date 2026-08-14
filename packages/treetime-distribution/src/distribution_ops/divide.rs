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
      let function = DistributionFunction::from_range_values(bounds, values)?;
      Ok(Distribution::Function(function))
    },
  }
}

fn division_support_intersection<Y: YAxisPolicy>(
  dividend_bounds: (f64, f64),
  divisor: &DistributionFunction<f64, Y>,
) -> SupportIntersection {
  let divisor_bounds = (
    if divisor.left_extrap() == BoundaryBehavior::Error {
      divisor.x_min()
    } else {
      dividend_bounds.0
    },
    if divisor.right_extrap() == BoundaryBehavior::Error {
      divisor.x_max()
    } else {
      dividend_bounds.1
    },
  );
  distribution_support_intersection(dividend_bounds, divisor_bounds)
}

/// The divisor's hard domain as division treats it, so [`guarded_empty_result`] reproduces exactly
/// the disjointness decided by [`division_support_intersection`].
///
/// A side bounds the quotient only where the divisor is undefined beyond the grid edge (`Error`).
/// Any other tail leaves the divisor defined-or-zero beyond (`Y::safe_divisor` absorbs the zero), so
/// the quotient continues under the dividend and that side does not separate the domains -- modelled
/// here as a soft (unbounded) boundary.
fn division_divisor_domain<Y: YAxisPolicy>(divisor: &DistributionFunction<f64, Y>) -> HardDomain {
  let bounding = |tail: BoundaryBehavior| {
    if tail == BoundaryBehavior::Error {
      BoundaryBehavior::Error
    } else {
      BoundaryBehavior::Constant
    }
  };
  (
    (divisor.x_min(), divisor.x_max()),
    (bounding(divisor.left_extrap()), bounding(divisor.right_extrap())),
  )
}
