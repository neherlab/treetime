use crate::Distribution;
use crate::distribution_core::formula::DistributionFormula;
use crate::distribution_core::function::DistributionFunction;
use crate::distribution_core::point::DistributionPoint;
use crate::distribution_core::range::DistributionRange;
use crate::distribution_ops::time_bounds::{
  SupportIntersection, distribution_support_intersection, distribution_support_n_points,
};
use crate::policy::YAxisPolicy;
use eyre::Report;
use ndarray::{Array1, Zip};
use treetime_grid::{ApproachLaw, BoundaryBehavior};
use treetime_utils::make_internal_error;

/// Grid size for discretizing Formula distributions that have no natural grid.
const FORMULA_GRID_SIZE: usize = 200;

pub fn distribution_multiplication<Y: YAxisPolicy>(
  a: &Distribution<Y>,
  b: &Distribution<Y>,
) -> Result<Distribution<Y>, Report> {
  match (a, b) {
    (Distribution::Empty, _) | (_, Distribution::Empty) => {
      Ok(Distribution::empty()) //
    },
    (Distribution::Point(a), Distribution::Point(b)) => {
      multiply_point_point::<Y>(a, b) //
    },
    (Distribution::Point(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Point(a)) => {
      multiply_point_function::<Y>(a, b) //
    },
    (Distribution::Point(a), Distribution::Range(b)) | (Distribution::Range(b), Distribution::Point(a)) => {
      multiply_point_range::<Y>(a, b) //
    },
    (Distribution::Range(a), Distribution::Range(b)) => {
      multiply_range_range::<Y>(a, b) //
    },
    (Distribution::Range(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Range(a)) => {
      multiply_range_function::<Y>(a, b) //
    },
    (Distribution::Function(a), Distribution::Function(b)) => {
      multiply_function_function::<Y>(a, b) //
    },
    (Distribution::Formula(a), Distribution::Formula(b)) => {
      multiply_formula_formula::<Y>(a, b) //
    },
    (Distribution::Formula(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Formula(a)) => {
      multiply_formula_function::<Y>(a, b) //
    },
    (Distribution::Formula(a), Distribution::Point(b)) | (Distribution::Point(b), Distribution::Formula(a)) => {
      multiply_formula_point::<Y>(a, b) //
    },
    (Distribution::Formula(a), Distribution::Range(b)) | (Distribution::Range(b), Distribution::Formula(a)) => {
      multiply_formula_range::<Y>(a, b) //
    },
  }
}

fn multiply_point_point<Y: YAxisPolicy>(
  a: &DistributionPoint<f64, Y>,
  b: &DistributionPoint<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  const EPS: f64 = 1e-9;
  if (a.t() - b.t()).abs() > EPS {
    return Ok(Distribution::empty());
  }
  let amplitude = Y::multiply(a.amplitude(), b.amplitude());
  Ok(Distribution::point(a.t(), amplitude))
}

fn multiply_point_range<Y: YAxisPolicy>(
  point: &DistributionPoint<f64, Y>,
  range: &DistributionRange<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  const EPS: f64 = 1e-9;
  let t = point.t();
  if t < range.start() - EPS || t > range.end() + EPS {
    return Ok(Distribution::empty());
  }
  let amplitude = Y::multiply(point.amplitude(), range.amplitude());
  Ok(Distribution::point(t, amplitude))
}

fn multiply_point_function<Y: YAxisPolicy>(
  point: &DistributionPoint<f64, Y>,
  func: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let t = point.t();
  let tail = if t < func.x_min() {
    Some(func.left_extrap())
  } else if t > func.x_max() {
    Some(func.right_extrap())
  } else {
    None
  };
  if let Some(tail) = tail {
    match tail {
      BoundaryBehavior::Constant => {},
      BoundaryBehavior::Error | BoundaryBehavior::Hard => return Ok(Distribution::empty()),
    }
  }
  let func_value = func.interp(t)?;
  let amplitude = Y::multiply(point.amplitude(), func_value);
  if !Y::is_defined(amplitude) {
    return Ok(Distribution::empty());
  }
  Ok(Distribution::point(t, amplitude))
}

fn multiply_range_range<Y: YAxisPolicy>(
  a: &DistributionRange<f64, Y>,
  b: &DistributionRange<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let overlap_start = a.start().max(b.start());
  let overlap_end = a.end().min(b.end());

  if overlap_start >= overlap_end {
    return Ok(Distribution::empty());
  }

  let amplitude = Y::multiply(a.amplitude(), b.amplitude());
  Ok(Distribution::range((overlap_start, overlap_end), amplitude))
}

fn multiply_range_function<Y: YAxisPolicy>(
  range: &DistributionRange<f64, Y>,
  func: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let a_bounds = (range.start(), range.end());
  let a_tails = (BoundaryBehavior::Error, BoundaryBehavior::Error);
  let b_bounds = (func.x_min(), func.x_max());
  let b_tails = (func.left_extrap(), func.right_extrap());
  match multiplication_support_intersection(a_bounds, a_tails, b_bounds, b_tails) {
    SupportIntersection::Disjoint => multiplication_empty_result(a_bounds, a_tails, b_bounds, b_tails),
    SupportIntersection::Point(t) => {
      let amplitude = Y::multiply(range.amplitude(), func.interp(t)?);
      Ok(Distribution::point(t, amplitude))
    },
    SupportIntersection::Interval(bounds) => {
      let n_points = distribution_support_n_points(bounds, func.dx())?;
      let grid = Array1::linspace(bounds.0, bounds.1, n_points);
      let values = func
        .interp_many(&grid)?
        .mapv(|value| Y::multiply(range.amplitude(), value));
      // A range has no interpolated tail (Error both sides), so composition keeps Error tails.
      let function = with_composed_tails(
        DistributionFunction::from_range_values(bounds, values)?,
        a_tails,
        b_tails,
      )?;
      Ok(Distribution::Function(function))
    },
  }
}

fn multiply_function_function<Y: YAxisPolicy>(
  a: &DistributionFunction<f64, Y>,
  b: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let a_bounds = (a.x_min(), a.x_max());
  let a_tails = (a.left_extrap(), a.right_extrap());
  let b_bounds = (b.x_min(), b.x_max());
  let b_tails = (b.left_extrap(), b.right_extrap());
  let a_approaches = (a.left_approach().copied(), a.right_approach().copied());
  let b_approaches = (b.left_approach().copied(), b.right_approach().copied());
  match multiplication_support_intersection(a_bounds, a_tails, b_bounds, b_tails) {
    SupportIntersection::Disjoint => multiplication_empty_result(a_bounds, a_tails, b_bounds, b_tails),
    SupportIntersection::Point(t) => Ok(Distribution::point(t, Y::multiply(a.interp(t)?, b.interp(t)?))),
    SupportIntersection::Interval(bounds) => {
      let n_points = distribution_support_n_points(bounds, a.dx().min(b.dx()))?;
      let grid = Array1::linspace(bounds.0, bounds.1, n_points);
      let values_a = a.interp_many(&grid)?;
      let values_b = b.interp_many(&grid)?;
      let values = Zip::from(&values_a)
        .and(&values_b)
        .map_collect(|&value_a, &value_b| Y::multiply(value_a, value_b));
      let function = with_composed_tails_and_approaches(
        DistributionFunction::from_range_values(bounds, values)?,
        a_tails,
        b_tails,
        a_approaches,
        b_approaches,
      )?;
      Ok(Distribution::Function(function))
    },
  }
}

fn multiply_formula_formula<Y: YAxisPolicy>(
  a: &DistributionFormula<Y>,
  b: &DistributionFormula<Y>,
) -> Result<Distribution<Y>, Report> {
  let overlap_min = a.t_min().max(b.t_min());
  let overlap_max = a.t_max().min(b.t_max());

  if overlap_min >= overlap_max {
    return Ok(Distribution::empty());
  }

  let n_points = FORMULA_GRID_SIZE;
  let values = (0..n_points)
    .map(|i| {
      let t = overlap_min + (overlap_max - overlap_min) * (i as f64 / (n_points - 1) as f64);
      let va = a.eval_single(t)?;
      let vb = b.eval_single(t)?;
      Ok(Y::multiply(va, vb))
    })
    .collect::<Result<Vec<f64>, Report>>()?;

  let distribution_fn = DistributionFunction::from_range_values((overlap_min, overlap_max), Array1::from_vec(values))?;
  Ok(Distribution::Function(distribution_fn))
}

fn multiply_formula_function<Y: YAxisPolicy>(
  a: &DistributionFormula<Y>,
  b: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let a_bounds = (a.t_min(), a.t_max());
  let a_tails = (BoundaryBehavior::Error, BoundaryBehavior::Error);
  let b_bounds = (b.x_min(), b.x_max());
  let b_tails = (b.left_extrap(), b.right_extrap());
  match multiplication_support_intersection(a_bounds, a_tails, b_bounds, b_tails) {
    SupportIntersection::Disjoint => multiplication_empty_result(a_bounds, a_tails, b_bounds, b_tails),
    SupportIntersection::Point(t) => Ok(Distribution::point(t, Y::multiply(a.eval_single(t)?, b.interp(t)?))),
    SupportIntersection::Interval(bounds) => {
      let n_points = distribution_support_n_points(bounds, b.dx())?;
      let grid = Array1::linspace(bounds.0, bounds.1, n_points);
      let formula_values = a.eval_many(&grid)?;
      let function_values = b.interp_many(&grid)?;
      let values = Zip::from(&formula_values)
        .and(&function_values)
        .map_collect(|&formula, &function| Y::multiply(formula, function));
      // A formula has no interpolated tail (Error both sides), so composition keeps Error tails.
      let function = with_composed_tails(
        DistributionFunction::from_range_values(bounds, values)?,
        a_tails,
        b_tails,
      )?;
      Ok(Distribution::Function(function))
    },
  }
}

fn multiply_formula_point<Y: YAxisPolicy>(
  a: &DistributionFormula<Y>,
  b: &DistributionPoint<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  const EPS: f64 = 1e-9;
  let t = b.t();

  if t < a.t_min() - EPS || t > a.t_max() + EPS {
    return Ok(Distribution::empty());
  }

  let va = a.eval_single(t)?;
  let amplitude = Y::multiply(va, b.amplitude());

  if !Y::is_defined(amplitude) {
    return Ok(Distribution::empty());
  }

  Ok(Distribution::point(t, amplitude))
}

fn multiply_formula_range<Y: YAxisPolicy>(
  a: &DistributionFormula<Y>,
  b: &DistributionRange<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let overlap_min = a.t_min().max(b.start());
  let overlap_max = a.t_max().min(b.end());

  if overlap_min >= overlap_max {
    return Ok(Distribution::empty());
  }

  let b_amplitude = b.amplitude();
  let n_points = FORMULA_GRID_SIZE;
  let values = (0..n_points)
    .map(|i| {
      let t = overlap_min + (overlap_max - overlap_min) * (i as f64 / (n_points - 1) as f64);
      let va = a.eval_single(t)?;
      Ok(Y::multiply(va, b_amplitude))
    })
    .collect::<Result<Vec<f64>, Report>>()?;

  let distribution_fn = DistributionFunction::from_range_values((overlap_min, overlap_max), Array1::from_vec(values))?;
  Ok(Distribution::Function(distribution_fn))
}

/// Attach the per-side result tails composed from the two operands' tails to a Function result.
fn with_composed_tails<Y: YAxisPolicy>(
  function: DistributionFunction<f64, Y>,
  a_tails: (BoundaryBehavior, BoundaryBehavior),
  b_tails: (BoundaryBehavior, BoundaryBehavior),
) -> Result<DistributionFunction<f64, Y>, Report> {
  function
    .with_left_extrap(compose_multiplication_tail(a_tails.0, b_tails.0))?
    .with_right_extrap(compose_multiplication_tail(a_tails.1, b_tails.1))
}

/// Attach composed tails and approach laws from both operands to a Function result.
fn with_composed_tails_and_approaches<Y: YAxisPolicy>(
  function: DistributionFunction<f64, Y>,
  a_tails: (BoundaryBehavior, BoundaryBehavior),
  b_tails: (BoundaryBehavior, BoundaryBehavior),
  a_approaches: (Option<ApproachLaw>, Option<ApproachLaw>),
  b_approaches: (Option<ApproachLaw>, Option<ApproachLaw>),
) -> Result<DistributionFunction<f64, Y>, Report> {
  let result = function
    .with_left_extrap(compose_multiplication_tail(a_tails.0, b_tails.0))?
    .with_right_extrap(compose_multiplication_tail(a_tails.1, b_tails.1))?;
  let left_approach = compose_approach_laws(a_approaches.0, b_approaches.0);
  let right_approach = compose_approach_laws(a_approaches.1, b_approaches.1);
  Ok(result.with_left_approach(left_approach).with_right_approach(right_approach))
}

/// Compose approach laws from two operands on the same side under multiplication.
///
/// Both present: exponents add, coefficients multiply.
/// One present: the operand without an approach law contributes its grid boundary
/// value as a constant factor (exponent 0), so the single law is returned as-is
/// (the grid boundary value is already accounted for in the product's grid values).
/// Neither present: no approach law on the result.
fn compose_approach_laws(a: Option<ApproachLaw>, b: Option<ApproachLaw>) -> Option<ApproachLaw> {
  match (a, b) {
    (Some(a), Some(b)) => Some(a.compose_multiply(&b)),
    (Some(law), None) | (None, Some(law)) => Some(law),
    (None, None) => None,
  }
}

/// Compose the per-side result tail for a product from the two operands' tails on that side.
///
/// Beyond a boundary the product is evaluated pointwise: if either operand is undefined there
/// the product is undefined (`Error`); otherwise if either operand is zero the product is zero
/// (`Hard`); only when both operands are flat non-zero constants is the product a flat constant
/// (`Constant`). This is the maximum over the precedence `Constant` < `Hard` < `Error`.
fn compose_multiplication_tail(a: BoundaryBehavior, b: BoundaryBehavior) -> BoundaryBehavior {
  match (a, b) {
    (BoundaryBehavior::Error, _) | (_, BoundaryBehavior::Error) => BoundaryBehavior::Error,
    (BoundaryBehavior::Hard, _) | (_, BoundaryBehavior::Hard) => BoundaryBehavior::Hard,
    (BoundaryBehavior::Constant, BoundaryBehavior::Constant) => BoundaryBehavior::Constant,
  }
}

/// Compute support intersection for multiplication, honoring operand tails.
///
/// Each grid boundary is either *hard* (the domain terminates: `Hard` is zero probability
/// beyond, `Error` is undefined beyond) or *soft* (the distribution continues past the grid
/// edge under a declared tail law). The product domain, resolved per side independently, is:
///
/// - both hard: the *tightest* (innermost) hard bound -- a hard bound is a fact about a
///   distribution, so the product can be non-zero only where both operands are;
/// - both soft: the *loosest* (outermost) soft bound -- either operand's tail law continues to
///   the other, so the product remains evaluable out to the farther edge;
/// - mixed: the hard bound dominates -- it terminates the domain regardless of the soft operand.
///
/// A soft boundary extends the operand's evaluable domain to the other operand's bound on the
/// same side; a hard boundary keeps its grid edge. Intersecting the two extended domains then
/// selects the tightest-hard / loosest-soft bound automatically. This is why the Ebola-scale
/// disjoint-grid product stays non-empty: the backward messages are soft on the left (the parent
/// could be arbitrarily older), so the left side takes the loosest soft bound instead of
/// collapsing to an empty intersection of the two finite grids.
fn multiplication_support_intersection(
  a_bounds: (f64, f64),
  a_tails: (BoundaryBehavior, BoundaryBehavior),
  b_bounds: (f64, f64),
  b_tails: (BoundaryBehavior, BoundaryBehavior),
) -> SupportIntersection {
  let a_eval = (
    if a_tails.0.is_soft() {
      a_bounds.0.min(b_bounds.0)
    } else {
      a_bounds.0
    },
    if a_tails.1.is_soft() {
      a_bounds.1.max(b_bounds.1)
    } else {
      a_bounds.1
    },
  );
  let b_eval = (
    if b_tails.0.is_soft() {
      b_bounds.0.min(a_bounds.0)
    } else {
      b_bounds.0
    },
    if b_tails.1.is_soft() {
      b_bounds.1.max(a_bounds.1)
    } else {
      b_bounds.1
    },
  );
  distribution_support_intersection(a_eval, b_eval)
}

/// Produce the `Empty` result of a Function-producing multiplication, checking that it is
/// legitimate.
///
/// A product is `Empty` only when the operands' *hard* domains are genuinely disjoint: a hard
/// boundary is a fact about a distribution (probability is zero, or evaluation undefined, beyond
/// it), so two distributions whose hard domains do not overlap cannot both be non-zero anywhere.
/// A soft boundary never separates domains -- the distribution continues past it under its tail
/// law -- so an `Empty` that disjoint hard domains do not explain is a numerical or logic
/// collapse. In the timetree passes such an `Empty` silently poisons every ancestor to the root
/// (the original motivating defect), so it is reported as an internal error rather than returned.
fn multiplication_empty_result<Y: YAxisPolicy>(
  a_bounds: (f64, f64),
  a_tails: (BoundaryBehavior, BoundaryBehavior),
  b_bounds: (f64, f64),
  b_tails: (BoundaryBehavior, BoundaryBehavior),
) -> Result<Distribution<Y>, Report> {
  if hard_domains_disjoint(a_bounds, a_tails, b_bounds, b_tails) {
    Ok(Distribution::empty())
  } else {
    make_internal_error!(
      "Multiplication produced an empty result, but the operands' hard domains overlap: \
       operand A grid [{}, {}] tails ({:?}, {:?}), operand B grid [{}, {}] tails ({:?}, {:?}). \
       An empty product must arise only from genuinely disjoint hard domains, never from numerical collapse",
      a_bounds.0,
      a_bounds.1,
      a_tails.0,
      a_tails.1,
      b_bounds.0,
      b_bounds.1,
      b_tails.0,
      b_tails.1
    )
  }
}

/// Whether the operands' *hard* domains are genuinely disjoint.
///
/// A soft boundary does not bound the domain (the distribution continues under its tail law), so
/// it is treated as unbounded on that side; only hard boundaries can separate two domains. The
/// hard domains are disjoint when the greatest hard lower bound exceeds the least hard upper
/// bound.
pub fn hard_domains_disjoint(
  a_bounds: (f64, f64),
  a_tails: (BoundaryBehavior, BoundaryBehavior),
  b_bounds: (f64, f64),
  b_tails: (BoundaryBehavior, BoundaryBehavior),
) -> bool {
  let a_lo = if a_tails.0.is_soft() {
    f64::NEG_INFINITY
  } else {
    a_bounds.0
  };
  let a_hi = if a_tails.1.is_soft() { f64::INFINITY } else { a_bounds.1 };
  let b_lo = if b_tails.0.is_soft() {
    f64::NEG_INFINITY
  } else {
    b_bounds.0
  };
  let b_hi = if b_tails.1.is_soft() { f64::INFINITY } else { b_bounds.1 };
  a_lo.max(b_lo) > a_hi.min(b_hi)
}
