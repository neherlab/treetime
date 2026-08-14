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
use treetime_utils::make_internal_error;

/// Grid size for discretizing Formula distributions that have no natural grid.
const FORMULA_GRID_SIZE: usize = 200;

pub fn distribution_multiplication<Y: YAxisPolicy>(
  a: &Distribution<Y>,
  b: &Distribution<Y>,
) -> Result<Distribution<Y>, Report> {
  match (a, b) {
    (Distribution::Empty, _) | (_, Distribution::Empty) => guarded_empty_result(
      "multiplication",
      distribution_hard_domain(a),
      distribution_hard_domain(b),
    ),
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
    return guarded_empty_result("multiplication", Some(point_hard_domain(a)), Some(point_hard_domain(b)));
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
    return guarded_empty_result(
      "multiplication",
      Some(point_hard_domain(point)),
      Some(range_hard_domain(range)),
    );
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
    // Beyond a soft boundary (`Constant` or `Linear`) the function continues under its tail law,
    // so the point picks up a finite value. Beyond a hard boundary (`Hard`) or an undeclared one
    // (`Error`) the function carries no probability there and the product is empty.
    if !tail.is_soft() {
      return guarded_empty_result(
        "multiplication",
        Some(point_hard_domain(point)),
        Some(function_hard_domain(func)),
      );
    }
  }
  let func_value = func.interp(t)?;
  let amplitude = Y::multiply(point.amplitude(), func_value);
  if !Y::is_defined(amplitude) {
    return guarded_empty_result(
      "multiplication",
      Some(point_hard_domain(point)),
      Some(function_hard_domain(func)),
    );
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
    return guarded_empty_result("multiplication", Some(range_hard_domain(a)), Some(range_hard_domain(b)));
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
    SupportIntersection::Disjoint => {
      guarded_empty_result("multiplication", Some((a_bounds, a_tails)), Some((b_bounds, b_tails)))
    },
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
  match multiplication_support_intersection(a_bounds, a_tails, b_bounds, b_tails) {
    SupportIntersection::Disjoint => {
      guarded_empty_result("multiplication", Some((a_bounds, a_tails)), Some((b_bounds, b_tails)))
    },
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
        a_tails,
        b_tails,
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
    return guarded_empty_result(
      "multiplication",
      Some(formula_hard_domain(a)),
      Some(formula_hard_domain(b)),
    );
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
    SupportIntersection::Disjoint => {
      guarded_empty_result("multiplication", Some((a_bounds, a_tails)), Some((b_bounds, b_tails)))
    },
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
    return guarded_empty_result(
      "multiplication",
      Some(formula_hard_domain(a)),
      Some(point_hard_domain(b)),
    );
  }

  let va = a.eval_single(t)?;
  let amplitude = Y::multiply(va, b.amplitude());

  if !Y::is_defined(amplitude) {
    return guarded_empty_result(
      "multiplication",
      Some(formula_hard_domain(a)),
      Some(point_hard_domain(b)),
    );
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
    return guarded_empty_result(
      "multiplication",
      Some(formula_hard_domain(a)),
      Some(range_hard_domain(b)),
    );
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
/// The product is pointwise, so the result class is the strongest of soft < `Hard` < `Error`:
/// `Error` if either operand is undefined beyond the edge; `Hard` if either operand ends the
/// domain (a hard bound restricts the product regardless of the other side); soft only when both
/// operands continue softly.
///
/// Fitted laws compose in closed form:
///
/// - Two `Hard` approach laws compose by adding all parameters (anchors, exponents, and slopes
///   add under multiplication in neg-log space). If only one operand
///   carries a law, the result carries none (`Hard(None)`): the `None` operand declares zero
///   density in the sub-grid gap `[t_hard, t_first)`, so the product is zero there and the present
///   law must not survive.
/// - Two `Linear` soft tails compose by adding their neg-log slopes (multiplication is addition in
///   neg-log space). A `Linear` tail times a flat `Constant` keeps the `Linear` slope, because a
///   flat tail contributes slope zero. `Linear(None)` times anything yields `Linear(None)`: the
///   product is soft but its slope is unknown until refit.
fn compose_multiplication_tail(a: BoundaryBehavior, b: BoundaryBehavior) -> BoundaryBehavior {
  match (a, b) {
    (BoundaryBehavior::Error, _) | (_, BoundaryBehavior::Error) => BoundaryBehavior::Error,
    (BoundaryBehavior::Hard(a_law), BoundaryBehavior::Hard(b_law)) => {
      let composed = match (a_law, b_law) {
        (Some(a), Some(b)) => Some(a.compose_multiply(&b)),
        // NOTE: a `None` operand contributes zero density in the sub-grid gap
        // `[t_hard, t_first)`, so the product is zero there. Never propagate the present
        // law -- that would fabricate density the `None` operand forbids.
        (Some(_) | None, None) | (None, Some(_)) => None,
      };
      BoundaryBehavior::Hard(composed)
    },
    (BoundaryBehavior::Hard(law), _) | (_, BoundaryBehavior::Hard(law)) => BoundaryBehavior::Hard(law),
    (BoundaryBehavior::Linear(a_law), BoundaryBehavior::Linear(b_law)) => {
      let composed = match (a_law, b_law) {
        (Some(a), Some(b)) => Some(a.compose_multiply(&b)),
        _ => None,
      };
      BoundaryBehavior::Linear(composed)
    },
    // A flat `Constant` contributes slope zero, so it leaves the `Linear` slope unchanged.
    (BoundaryBehavior::Linear(law), BoundaryBehavior::Constant)
    | (BoundaryBehavior::Constant, BoundaryBehavior::Linear(law)) => BoundaryBehavior::Linear(law),
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

/// A distribution's hard support domain: grid bounds `(lo, hi)` and the boundary behavior on each
/// side. `None` (see [`distribution_hard_domain`]) denotes an empty support -- the empty set, which
/// is disjoint from every domain.
pub type HardDomain = ((f64, f64), (BoundaryBehavior, BoundaryBehavior));

/// Hard domain of a point mass: the single point `t`, hard on both sides.
pub fn point_hard_domain<Y: YAxisPolicy>(p: &DistributionPoint<f64, Y>) -> HardDomain {
  ((p.t(), p.t()), (BoundaryBehavior::Error, BoundaryBehavior::Error))
}

/// Hard domain of a range: its finite support, hard on both sides.
pub fn range_hard_domain<Y: YAxisPolicy>(r: &DistributionRange<f64, Y>) -> HardDomain {
  ((r.start(), r.end()), (BoundaryBehavior::Error, BoundaryBehavior::Error))
}

/// Hard domain of a function: its grid bounds and the declared per-side tail behavior.
pub fn function_hard_domain<Y: YAxisPolicy>(f: &DistributionFunction<f64, Y>) -> HardDomain {
  ((f.x_min(), f.x_max()), (f.left_extrap(), f.right_extrap()))
}

/// Hard domain of a formula: its evaluable range, hard on both sides.
pub fn formula_hard_domain<Y: YAxisPolicy>(f: &DistributionFormula<Y>) -> HardDomain {
  (
    (f.t_min(), f.t_max()),
    (BoundaryBehavior::Error, BoundaryBehavior::Error),
  )
}

/// Hard support domain of a distribution, or `None` when its support is empty.
///
/// An empty support has no bounds and is disjoint from every other domain, so `None` marks the
/// legitimate-empty case for [`guarded_empty_result`].
pub fn distribution_hard_domain<Y: YAxisPolicy>(d: &Distribution<Y>) -> Option<HardDomain> {
  match d {
    Distribution::Empty => None,
    Distribution::Point(p) => Some(point_hard_domain(p)),
    Distribution::Range(r) => Some(range_hard_domain(r)),
    Distribution::Formula(f) => Some(formula_hard_domain(f)),
    Distribution::Function(f) => (!f.is_empty()).then(|| function_hard_domain(f)),
  }
}

/// Produce a guarded `Empty` result, verifying that the emptiness is legitimate.
///
/// `Distribution::Empty` must be reachable only from genuinely disjoint hard domains, or from an
/// operand whose support is already empty (`None`). A hard boundary is a fact about a distribution
/// (probability is zero, or evaluation undefined, beyond it), so two distributions whose hard
/// domains do not overlap cannot both carry mass anywhere; a soft boundary never separates domains,
/// because the distribution continues past it under its tail law. An `Empty` that disjoint hard
/// domains do not explain is a numerical or logic collapse. In the timetree passes such an `Empty`
/// silently poisons every ancestor to the root (the original motivating defect), so it is reported
/// as an internal error rather than returned. `operation` names the operation for the diagnostic.
pub fn guarded_empty_result<Y: YAxisPolicy>(
  operation: &str,
  a: Option<HardDomain>,
  b: Option<HardDomain>,
) -> Result<Distribution<Y>, Report> {
  let disjoint = match (a, b) {
    // An operand with empty support contributes the empty set, disjoint from every domain.
    (None, _) | (_, None) => true,
    (Some((a_bounds, a_tails)), Some((b_bounds, b_tails))) => {
      hard_domains_disjoint(a_bounds, a_tails, b_bounds, b_tails)
    },
  };
  if disjoint {
    Ok(Distribution::empty())
  } else {
    make_internal_error!(
      "{operation} produced an empty result, but the operands' hard domains overlap: \
       operand A {a:?}, operand B {b:?}. An empty result must arise only from genuinely disjoint \
       hard domains, never from numerical collapse"
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
