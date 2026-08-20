use crate::Distribution;
use crate::distribution_core::formula::DistributionFormula;
use crate::distribution_core::function::DistributionFunction;
use crate::distribution_core::point::DistributionPoint;
use crate::distribution_core::range::DistributionRange;
use crate::distribution_ops::time_bounds::{SupportIntersection, distribution_support_n_points};
use crate::policy::YAxisPolicy;
use eyre::Report;
use itertools::izip;
use ndarray::{Array1, Zip};
use ordered_float::OrderedFloat;
use treetime_grid::{BoundaryBehavior, Side};
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
  match multiplication_support_intersection(&[(a_bounds, a_tails), (b_bounds, b_tails)]) {
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
  // Pairwise multiplication is the two-operand case of the N-ary co-location: co-locate both operands
  // on one common grid over the support intersection and reduce once. Routing both arities through the
  // single primitive keeps the grid construction identical, so `a * b` equals `product([a, b])`.
  multiply_functions(&[a, b])
}

/// Product of N gridded densities co-located on one common grid: the shared core of both
/// [`multiply_function_function`] and [`distribution_product`](crate::distribution_ops::product::distribution_product).
///
/// A product of independent factors is pointwise -- a sum of neg-log ordinates -- so it is exact,
/// commutative, and independent of factor order. All operands are co-located on one working grid over
/// the support intersection and interpolated once each (regardless of fan-out), then reduced by a
/// single elementwise product.
///
/// The grid follows the multiplication support contract (`kb/decisions/distribution-tails-and-arithmetic.md`):
/// per side the tightest (innermost) hard bound when any operand terminates the domain there,
/// otherwise the loosest (outermost) soft bound; spacing from the finest operand; both analytical
/// endpoints included via `Array1::linspace`, so a hard bound lands exactly on a grid node. A disjoint
/// intersection is a legitimate `Empty` (only hard sides can separate domains); endpoint-only contact
/// collapses to a `Point`.
///
/// Result tails compose in closed form by folding [`compose_multiplication_tail`] over the operands:
/// soft slopes add, a hard bound dominates, `Error` dominates all. No tail is re-fit from the summed
/// grid.
pub fn multiply_functions<Y: YAxisPolicy>(
  functions: &[&DistributionFunction<f64, Y>],
) -> Result<Distribution<Y>, Report> {
  // Float multiplication is commutative but not associative, so the reduction is folded in a canonical
  // operand order. This makes the product genuinely independent of caller/child order (bit-identical
  // across permutations of the same operands), not merely equal to rounding tolerance.
  let ordered = canonical_operand_order(functions);
  let (&first, rest) = ordered
    .split_first()
    .expect("multiply_functions requires at least one operand");

  let domains: Vec<HardDomain> = ordered.iter().copied().map(function_hard_domain).collect();
  match multiplication_support_intersection(&domains) {
    // Genuinely disjoint hard domains carry no common support, so their product is legitimately empty.
    // A soft side never separates domains; only hard sides can, so this branch is the guarded-empty
    // case by construction.
    SupportIntersection::Disjoint => Ok(Distribution::empty()),
    // Supports meet at one point: the product collapses to a point mass whose amplitude is the
    // pointwise product of the operands sampled there (v0 converts the surviving knot to a delta).
    SupportIntersection::Point(t) => {
      let mut amplitude = first.interp(t)?;
      for f in rest {
        amplitude = Y::multiply(amplitude, f.interp(t)?);
      }
      Ok(Distribution::point(t, amplitude))
    },
    SupportIntersection::Interval(bounds) => {
      let dx = ordered
        .iter()
        .map(|f| f.dx())
        .reduce(f64::min)
        .expect("multiply_functions requires at least one operand");
      let n_points = distribution_support_n_points(bounds, dx)?;
      let grid = Array1::linspace(bounds.0, bounds.1, n_points);

      // Resample every operand onto the shared grid; each lands on the same grid, so the fold is a
      // plain elementwise product (a sum of ordinates in neg-log).
      let mut values = first.interp_many(&grid)?;
      for f in rest {
        let other = f.interp_many(&grid)?;
        Zip::from(&mut values)
          .and(&other)
          .for_each(|value, &o| *value = Y::multiply(*value, o));
      }

      let left_tail = compose_product_tail(ordered.iter().map(|f| f.left_extrap()))?;
      let right_tail = compose_product_tail(ordered.iter().map(|f| f.right_extrap()))?;
      let function = DistributionFunction::from_range_values(bounds, values)?
        .with_left_extrap(left_tail)?
        .with_right_extrap(right_tail)?;
      Ok(Distribution::Function(function))
    },
  }
}

/// Canonical, order-independent operand sequence for the reduction. Sorted by grid extent then
/// spacing, so any permutation of the same operands folds in the same order and yields a bit-identical
/// product despite floating-point non-associativity.
fn canonical_operand_order<'a, Y: YAxisPolicy>(
  functions: &[&'a DistributionFunction<f64, Y>],
) -> Vec<&'a DistributionFunction<f64, Y>> {
  let mut ordered = functions.to_vec();
  ordered.sort_by_key(|f| (OrderedFloat(f.x_min()), OrderedFloat(f.x_max()), OrderedFloat(f.dx())));
  ordered
}

/// Resolve the multiplication support intersection over the operands' hard domains: `Disjoint`, a
/// single-point `Point`, or an `Interval`.
///
/// Each grid boundary is either *hard* (the domain terminates: `Hard` is zero probability beyond,
/// `Error` is undefined beyond) or *soft* (the distribution continues past the edge under a declared
/// tail law). The product domain is resolved per side independently:
///
/// - any hard operand on the side: the tightest (innermost) hard bound -- a hard bound is a fact
///   about a distribution, so the product can be non-zero only where every hard operand is;
/// - otherwise: the loosest (outermost) soft bound -- each operand's tail law continues to the
///   others, so the product stays evaluable out to the farthest edge.
///
/// Disjointness is decided from hard sides only; a soft side continues under its tail law and never
/// separates domains. This is why the Ebola-scale disjoint-grid product stays non-empty: the backward
/// messages are soft on the left, so the left side takes the loosest soft bound instead of collapsing
/// to an empty intersection of the two finite grids. Endpoint contact uses exact comparison; a
/// tolerance would enlarge the intersection.
#[allow(clippy::float_cmp)] // Endpoint contact requires exact bound equality; a tolerance would enlarge the intersection.
fn multiplication_support_intersection(domains: &[HardDomain]) -> SupportIntersection {
  // Resolve each side independently: hard operands take the tightest (innermost) bound and soft
  // operands the loosest (outermost). `None` means no operand of that class bounds the side.
  let (hard_lo, soft_lo) = side_bounds(domains, Side::Left);
  let (hard_hi, soft_hi) = side_bounds(domains, Side::Right);

  // Disjointness is a fact about hard sides only; a soft side is treated as unbounded on that side.
  // `izip!` over the two `Option`s yields a pair only when both hard bounds are present.
  if izip!(hard_lo, hard_hi).any(|(lo, hi)| lo > hi) {
    return SupportIntersection::Disjoint;
  }

  // A present hard bound overrides the soft bound on its side. With no operands each side keeps its
  // empty sentinel (lo = +inf, hi = -inf), which the final `else` reports as `Disjoint`.
  let lo = hard_lo.or(soft_lo).unwrap_or(f64::INFINITY);
  let hi = hard_hi.or(soft_hi).unwrap_or(f64::NEG_INFINITY);
  if lo == hi {
    SupportIntersection::Point(lo)
  } else if lo < hi {
    SupportIntersection::Interval((lo, hi))
  } else {
    // lo > hi while the hard domains overlap is a numerical collapse, not a real intersection; treat
    // it as disjoint so the caller yields a guarded empty rather than an inverted grid. Unreachable
    // for well-formed operands (lo > hi implies both sides hard, which the disjoint check caught).
    SupportIntersection::Disjoint
  }
}

/// Reduce the operand hard domains on one `side` into its hard and soft bound components. Hard
/// operands are combined into the tightest (innermost) bound and soft operands into the loosest
/// (outermost); either component is `None` when no operand of that class bounds the side. The
/// innermost and outermost reductions flip between the lower (`Left`) and upper (`Right`) side.
fn side_bounds(domains: &[HardDomain], side: Side) -> (Option<f64>, Option<f64>) {
  let (inner, outer): (fn(f64, f64) -> f64, fn(f64, f64) -> f64) = match side {
    Side::Left => (f64::max, f64::min),
    Side::Right => (f64::min, f64::max),
  };
  let extrap = |d: &HardDomain| match side {
    Side::Left => d.1.0,
    Side::Right => d.1.1,
  };
  let bound = |d: &HardDomain| match side {
    Side::Left => d.0.0,
    Side::Right => d.0.1,
  };
  let hard = domains
    .iter()
    .filter(|&d| !extrap(d).is_soft())
    .map(&bound)
    .reduce(inner);
  let soft = domains
    .iter()
    .filter(|&d| extrap(d).is_soft())
    .map(&bound)
    .reduce(outer);
  (hard, soft)
}

/// Compose the per-side result tail of an N-ary product by folding the pairwise composition
/// ([`compose_multiplication_tail`]): soft slopes add, a hard bound dominates, `Error` dominates all.
/// Closed-form, with no re-fit from the summed grid.
fn compose_product_tail(tails: impl Iterator<Item = BoundaryBehavior>) -> Result<BoundaryBehavior, Report> {
  let mut composed: Option<BoundaryBehavior> = None;
  for tail in tails {
    composed = Some(match composed {
      None => tail,
      Some(current) => compose_multiplication_tail(current, tail)?,
    });
  }
  Ok(composed.expect("compose_product_tail requires at least one operand"))
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
  match multiplication_support_intersection(&[(a_bounds, a_tails), (b_bounds, b_tails)]) {
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
    .with_left_extrap(compose_multiplication_tail(a_tails.0, b_tails.0)?)?
    .with_right_extrap(compose_multiplication_tail(a_tails.1, b_tails.1)?)
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
/// - A nullary `Hard` operand declares zero density in the sub-grid gap `[t_hard, t_first)`, so the
///   product is zero there: it absorbs a `HardApproach` law rather than propagating it, which would
///   fabricate density the `Hard` operand forbids.
/// - Two `Linear` soft tails compose by adding their neg-log slopes (multiplication is addition in
///   neg-log space).
///
/// Two `HardApproach` laws are not composed. At a shared boundary their product is the single power
/// law `b_1 + b_2`, but at distinct boundaries the tighter bound dominates while the other operand
/// contributes a smooth, locally varying factor there, which no single-exponent
/// [`HardApproachLaw`](treetime_grid::HardApproachLaw) at the tighter bound can represent. No
/// production path multiplies two hard-approach tails (message hard sides are the nullary `Hard`, and
/// only the branch-length factor produces `HardApproach`), so reaching this is a loud internal error
/// rather than a silent lossy composition.
pub(super) fn compose_multiplication_tail(
  a: BoundaryBehavior,
  b: BoundaryBehavior,
) -> Result<BoundaryBehavior, Report> {
  match (a, b) {
    (BoundaryBehavior::Error, _) | (_, BoundaryBehavior::Error) => Ok(BoundaryBehavior::Error),
    // Both sides hard: a nullary `Hard` zeroes the sub-grid gap and so absorbs a `HardApproach` law.
    (BoundaryBehavior::Hard, BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_))
    | (BoundaryBehavior::HardApproach(_), BoundaryBehavior::Hard) => Ok(BoundaryBehavior::Hard),
    (BoundaryBehavior::HardApproach(_), BoundaryBehavior::HardApproach(_)) => make_internal_error!(
      "Cannot multiply two HardApproach tails: their product is not representable by a \
       single-parameter hard-approach law, and this composition is unreachable in the inference pipeline"
    ),
    // A hard bound restricts the product regardless of the soft other side, so keep the hard side.
    (hard @ (BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_)), _) => Ok(hard),
    (_, hard @ (BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_))) => Ok(hard),
    (BoundaryBehavior::Linear(a_law), BoundaryBehavior::Linear(b_law)) => {
      Ok(BoundaryBehavior::Linear(a_law.compose_multiply(&b_law)))
    },
  }
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
