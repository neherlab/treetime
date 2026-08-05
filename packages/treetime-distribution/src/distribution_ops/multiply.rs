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
use treetime_grid::BoundaryBehavior;

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
      BoundaryBehavior::Error | BoundaryBehavior::Zero => return Ok(Distribution::empty()),
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
  match multiplication_support_intersection(
    (range.start(), range.end()),
    (BoundaryBehavior::Error, BoundaryBehavior::Error),
    (func.x_min(), func.x_max()),
    (func.left_extrap(), func.right_extrap()),
  ) {
    SupportIntersection::Disjoint => Ok(Distribution::empty()),
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
        (BoundaryBehavior::Error, BoundaryBehavior::Error),
        (func.left_extrap(), func.right_extrap()),
      )?;
      Ok(Distribution::Function(function))
    },
  }
}

fn multiply_function_function<Y: YAxisPolicy>(
  a: &DistributionFunction<f64, Y>,
  b: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  match multiplication_support_intersection(
    (a.x_min(), a.x_max()),
    (a.left_extrap(), a.right_extrap()),
    (b.x_min(), b.x_max()),
    (b.left_extrap(), b.right_extrap()),
  ) {
    SupportIntersection::Disjoint => Ok(Distribution::empty()),
    SupportIntersection::Point(t) => Ok(Distribution::point(t, Y::multiply(a.interp(t)?, b.interp(t)?))),
    SupportIntersection::Interval(bounds) => {
      let n_points = distribution_support_n_points(bounds, a.dx().min(b.dx()))?;
      let grid = Array1::linspace(bounds.0, bounds.1, n_points);
      let values_a = a.interp_many(&grid)?;
      let values_b = b.interp_many(&grid)?;
      let values = Zip::from(&values_a)
        .and(&values_b)
        .map_collect(|&value_a, &value_b| Y::multiply(value_a, value_b));
      let function = with_composed_tails(
        DistributionFunction::from_range_values(bounds, values)?,
        (a.left_extrap(), a.right_extrap()),
        (b.left_extrap(), b.right_extrap()),
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
  match multiplication_support_intersection(
    (a.t_min(), a.t_max()),
    (BoundaryBehavior::Error, BoundaryBehavior::Error),
    (b.x_min(), b.x_max()),
    (b.left_extrap(), b.right_extrap()),
  ) {
    SupportIntersection::Disjoint => Ok(Distribution::empty()),
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
        (BoundaryBehavior::Error, BoundaryBehavior::Error),
        (b.left_extrap(), b.right_extrap()),
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

/// Compose the per-side result tail for a product from the two operands' tails on that side.
///
/// Beyond a boundary the product is evaluated pointwise: if either operand is undefined there
/// the product is undefined (`Error`); otherwise if either operand is zero the product is zero
/// (`Zero`); only when both operands are flat non-zero constants is the product a flat constant
/// (`Constant`). This is the maximum over the precedence `Constant` < `Zero` < `Error`.
fn compose_multiplication_tail(a: BoundaryBehavior, b: BoundaryBehavior) -> BoundaryBehavior {
  match (a, b) {
    (BoundaryBehavior::Error, _) | (_, BoundaryBehavior::Error) => BoundaryBehavior::Error,
    (BoundaryBehavior::Zero, _) | (_, BoundaryBehavior::Zero) => BoundaryBehavior::Zero,
    (BoundaryBehavior::Constant, BoundaryBehavior::Constant) => BoundaryBehavior::Constant,
  }
}

/// Compute support intersection for multiplication, honoring operand tails.
///
/// Each grid boundary is either *hard* (the domain terminates: `Zero` is zero probability
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
