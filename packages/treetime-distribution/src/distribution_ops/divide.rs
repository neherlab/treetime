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
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, Side};
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
  // Beyond a hard divisor edge the divisor is zero (`Hard`/`HardApproach`, undefined quotient) or
  // itself undefined (`Error`), so the point carries no dividable mass there. Only a soft divisor
  // tail is evaluable and dividable beyond the grid. This mirrors `fn multiply_point_function()`.
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
) -> Result<Distribution<Y>, Report> {
  let range_domain = range_hard_domain(range);
  let divisor_domain = function_hard_domain(divisor);
  match multiplication_support_intersection(&[range_domain, divisor_domain]) {
    SupportIntersection::Disjoint => guarded_empty_result("division", Some(range_domain), Some(divisor_domain)),
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
      // A range is hard on both sides, so each result side inherits the range's `Hard` bound (the
      // quotient is zero outside the box), or `Error` where a divisor hard edge binds.
      let function = apply_division_tail(function, Side::Left, range_domain, divisor_domain)?;
      let function = apply_division_tail(function, Side::Right, range_domain, divisor_domain)?;
      Ok(Distribution::Function(function))
    },
  }
}

fn divide_function_by_function<Y: YAxisPolicy>(
  dividend: &DistributionFunction<f64, Y>,
  divisor: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  // Division shares multiplication's support rule: per side the innermost hard bound when any operand
  // is hard there, else the outermost soft bound. Both operands are evaluated through their own tail
  // laws on the common grid and their neg-log ordinates subtracted (division is subtraction in
  // neg-log). A soft divisor edge extends the quotient: the divisor's decaying tail is sampled as
  // bulk, and the dividend's own tail carries the cavity out to the union edge. The cavity
  // dividend is a product that contains the divisor as a factor, so the quotient's tail decays rather
  // than exploding. See `kb/decisions/distribution-tails-and-arithmetic.md` (Division).
  let dividend_domain = function_hard_domain(dividend);
  let divisor_domain = function_hard_domain(divisor);
  match multiplication_support_intersection(&[dividend_domain, divisor_domain]) {
    SupportIntersection::Disjoint => guarded_empty_result("division", Some(dividend_domain), Some(divisor_domain)),
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
      let function = apply_division_tail(function, Side::Left, dividend_domain, divisor_domain)?;
      let function = apply_division_tail(function, Side::Right, dividend_domain, divisor_domain)?;
      Ok(Distribution::Function(function))
    },
  }
}

/// Attach the per-side result tail of a Function quotient, refitting a soft side from the built grid.
fn apply_division_tail<Y: YAxisPolicy>(
  function: DistributionFunction<f64, Y>,
  side: Side,
  dividend: HardDomain,
  divisor: HardDomain,
) -> Result<DistributionFunction<f64, Y>, Report> {
  match division_side_tail(side, dividend, divisor) {
    // Both operands soft: refit the quotient's own decaying `Linear` tail from the combined grid,
    // matching the Python prototype `combine()`. Exact slope subtraction holds only while every
    // operand's grid reaches the same edge, which the cavity's independently regridded operands do
    // not, so the tail is refit rather than composed in closed form. `SoftTailLaw::fit` clamps a
    // wrong-sign fit to a flat tail, which is the guard against a tail that manufactures mass outward.
    DivisionSideTail::Refit => function.fit_soft_tail(side, DEFAULT_TAIL_FIT_POINTS),
    DivisionSideTail::Fixed(tail) => match side {
      Side::Left => function.with_left_extrap(tail),
      Side::Right => function.with_right_extrap(tail),
    },
  }
}

/// Per-side result tail of a Function quotient, deferred where a soft side must be refit from the
/// built grid.
///
/// The grid rule is multiplication's (intersect hard, union soft), so the intersection edge on a side
/// is the innermost hard bound when any operand is hard there, else the outermost soft bound.
/// Division's tail differs from multiplication's only in the hard cases:
///
/// - **both soft** -> the divisor's soft tail was sampled as bulk out to the union edge; beyond it the
///   quotient continues, so refit a decaying `Linear` law from the quotient grid.
/// - **only the dividend hard** -> the dividend binds; beyond its edge the dividend is zero
///   (`Hard`/`HardApproach`) so the quotient is zero, or undefined (`Error`) so the quotient is
///   undefined. Inherit the dividend's own hard tail.
/// - **only the divisor hard** -> the divisor binds; beyond its edge the divisor is zero (a `f/0`
///   spike) or undefined, so the quotient is `Error`.
/// - **both hard** -> the innermost bound binds; a tie is a `0/0` cavity edge the dividend owns (the
///   cavity dividend is a product containing the divisor as a factor), so the dividend edge inherits
///   the dividend tail and the divisor edge is `Error`.
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

/// The [`BoundaryBehavior`] a hard domain declares on one `side`.
fn side_tail(side: Side, domain: HardDomain) -> BoundaryBehavior {
  match side {
    Side::Left => domain.1.0,
    Side::Right => domain.1.1,
  }
}

/// The grid bound a hard domain carries on one `side`.
fn side_bound(side: Side, domain: HardDomain) -> f64 {
  match side {
    Side::Left => domain.0.0,
    Side::Right => domain.0.1,
  }
}

/// Outcome of resolving one quotient side: refit a soft tail, or set a fixed hard one.
enum DivisionSideTail {
  /// Both operands are soft on this side: refit a decaying `Linear` law from the quotient grid.
  Refit,
  /// The side is hard: the dividend's own hard tail when the dividend binds the edge, or `Error`
  /// when a divisor hard edge binds (the quotient is undefined or a spike past it).
  Fixed(BoundaryBehavior),
}
