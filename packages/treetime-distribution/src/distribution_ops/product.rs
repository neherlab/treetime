use crate::Distribution;
use crate::distribution_core::function::DistributionFunction;
use crate::distribution_ops::multiply::{compose_multiplication_tail, distribution_multiplication};
use crate::policy::YAxisPolicy;
use eyre::Report;
use ndarray::Zip;
use std::cmp::Ordering;
use treetime_grid::BoundaryBehavior;
use treetime_utils::make_internal_error;

/// Product of N distribution factors: the N-ary generalization of [`distribution_multiplication`].
///
/// A product of independent factors is pointwise -- a sum of neg-log ordinates -- so it is exact,
/// associative, and independent of factor order. The gridded `Function` factors are co-located on
/// one common grid and reduced once ([`multiply_functions`]), so every function is interpolated at
/// most once regardless of fan-out; no accumulator is re-resampled per factor. The remaining exact
/// factors (`Point`, `Range`, `Formula`) need no grid: their product is formed with the pairwise
/// [`distribution_multiplication`] and multiplied into the function product in one final step.
///
/// An `Empty` factor makes the product `Empty` (an empty operand is disjoint from every domain).
pub fn distribution_product<Y: YAxisPolicy>(factors: &[&Distribution<Y>]) -> Result<Distribution<Y>, Report> {
  let mut functions = Vec::new();
  let mut others = Vec::new();
  for &factor in factors {
    match factor {
      Distribution::Empty => return Ok(Distribution::empty()),
      Distribution::Function(function) => functions.push(function),
      Distribution::Point(_) | Distribution::Range(_) | Distribution::Formula(_) => others.push(factor),
    }
  }

  let function_product = if functions.is_empty() {
    None
  } else {
    Some(multiply_functions(&functions)?)
  };

  let mut other_product: Option<Distribution<Y>> = None;
  for factor in others {
    other_product = Some(match other_product {
      None => factor.clone(),
      Some(current) => distribution_multiplication(&current, factor)?,
    });
  }

  match (function_product, other_product) {
    (Some(functions), Some(others)) => distribution_multiplication(&functions, &others),
    (Some(single), None) | (None, Some(single)) => Ok(single),
    (None, None) => make_internal_error!("distribution_product requires at least one factor"),
  }
}

/// Multiply N gridded densities on one common grid: the N-ary form of the function product.
///
/// The product of gridded densities is pointwise (a sum of neg-log ordinates), so all operands are
/// co-located on one shared working grid and reduced once. The grid follows the multiplication
/// support rule per side: the tightest hard bound when any operand terminates the domain there,
/// otherwise the loosest soft bound; spacing from the finest operand. Each operand is resampled
/// exactly once onto that grid (tail-preserving), so no operand passes through more than one
/// interpolation regardless of fan-out, and the elementwise sum is independent of operand order.
///
/// The working grid keeps the finest operand's `dx` spacing (via
/// [`DistributionFunction::resample_range_dx`]); this is the grid the timetree backward child fold
/// established, distinct from the pairwise `multiply_function_function`'s inclusive `linspace`
/// intersection grid. The two agree to within a resolution step and both stay at least as fine as
/// their sharpest operand.
///
/// Result tails compose in closed form from the operand tails ([`compose_multiplication_tail`], the
/// N-ary generalization of `with_composed_tails`): soft `Linear` slopes add, a hard bound dominates,
/// `Error` dominates all. No tail is re-fit from the summed grid.
fn multiply_functions<Y: YAxisPolicy>(functions: &[&DistributionFunction<f64, Y>]) -> Result<Distribution<Y>, Report> {
  // Genuinely disjoint hard domains carry no common support, so their product is legitimately empty
  // (the guarded-empty invariant). A soft side never separates domains; only hard sides can.
  if functions_hard_domains_disjoint(functions) {
    return Ok(Distribution::empty());
  }

  let left = product_lower_bound(functions.iter().map(|f| (f.x_min(), f.left_extrap())));
  let right = product_upper_bound(functions.iter().map(|f| (f.x_max(), f.right_extrap())));
  let dx = functions
    .iter()
    .map(|f| f.dx())
    .reduce(f64::min)
    .expect("multiply_functions requires at least one operand");

  if !matches!(left.partial_cmp(&right), Some(Ordering::Less)) {
    // Overlapping hard domains (checked above) always leave left < right, so a non-`Less` compare --
    // equal, reversed, or `NaN` -- is a numerical collapse. Guarding it keeps an empty product from
    // silently poisoning every ancestor node.
    return make_internal_error!("distribution product produced an empty working grid [{left}, {right}]");
  }

  // Resample every operand onto the shared [left, right] / dx grid; each lands on the same grid, so
  // the fold is a plain elementwise product (a sum of ordinates in neg-log).
  let resampled = functions
    .iter()
    .map(|f| f.resample_range_dx((left, right), dx))
    .collect::<Result<Vec<_>, Report>>()?;
  let (first, rest) = resampled
    .split_first()
    .expect("multiply_functions requires at least one operand");

  let mut values = first.y().clone();
  for function in rest {
    Zip::from(&mut values)
      .and(function.y())
      .for_each(|value, &other| *value = Y::multiply(*value, other));
  }

  // Rebuild on the resampled grid's own extent (`resample_range_dx` ends at `left + (n - 1) * dx`,
  // a rounding step short of `right`), so the result keeps the exact grid the operands landed on.
  let left_tail = compose_product_tail(functions.iter().map(|f| f.left_extrap()))?;
  let right_tail = compose_product_tail(functions.iter().map(|f| f.right_extrap()))?;
  let function = DistributionFunction::from_range_values((first.x_min(), first.x_max()), values)?
    .with_left_extrap(left_tail)?
    .with_right_extrap(right_tail)?;
  Ok(Distribution::Function(function))
}

/// Lower (left) working-grid bound for a product: the tightest (innermost) hard bound when any
/// operand terminates its domain on the left, otherwise the loosest (outermost) soft bound. The
/// N-ary generalization of the left column of the pairwise `multiplication_support_intersection`.
fn product_lower_bound(operands: impl Iterator<Item = (f64, BoundaryBehavior)>) -> f64 {
  let mut hard: Option<f64> = None;
  let mut soft: Option<f64> = None;
  for (x_min, tail) in operands {
    if tail.is_soft() {
      soft = Some(soft.map_or(x_min, |s| s.min(x_min)));
    } else {
      hard = Some(hard.map_or(x_min, |h| h.max(x_min)));
    }
  }
  hard.or(soft).expect("at least one operand")
}

/// Upper (right) working-grid bound for a product: the tightest hard bound when any operand
/// terminates its domain on the right, otherwise the loosest soft bound.
fn product_upper_bound(operands: impl Iterator<Item = (f64, BoundaryBehavior)>) -> f64 {
  let mut hard: Option<f64> = None;
  let mut soft: Option<f64> = None;
  for (x_max, tail) in operands {
    if tail.is_soft() {
      soft = Some(soft.map_or(x_max, |s| s.max(x_max)));
    } else {
      hard = Some(hard.map_or(x_max, |h| h.min(x_max)));
    }
  }
  hard.or(soft).expect("at least one operand")
}

/// Whether the operands' hard domains are genuinely disjoint (no common support). A soft side does
/// not bound the domain -- the density continues under its tail law -- so it is treated as
/// unbounded; only hard sides can separate domains. The N-ary generalization of the pairwise
/// `hard_domains_disjoint`.
fn functions_hard_domains_disjoint<Y: YAxisPolicy>(functions: &[&DistributionFunction<f64, Y>]) -> bool {
  let mut hard_lo = f64::NEG_INFINITY;
  let mut hard_hi = f64::INFINITY;
  for f in functions {
    if !f.left_extrap().is_soft() {
      hard_lo = hard_lo.max(f.x_min());
    }
    if !f.right_extrap().is_soft() {
      hard_hi = hard_hi.min(f.x_max());
    }
  }
  hard_lo > hard_hi
}

/// Compose the per-side result tail of an N-ary product by folding the pairwise composition
/// ([`compose_multiplication_tail`]): soft `Linear` slopes add, a hard bound dominates, `Error`
/// dominates all. Closed-form, with no re-fit from the summed grid.
fn compose_product_tail(tails: impl Iterator<Item = BoundaryBehavior>) -> Result<BoundaryBehavior, Report> {
  let mut composed: Option<BoundaryBehavior> = None;
  for tail in tails {
    composed = Some(match composed {
      None => tail,
      Some(current) => compose_multiplication_tail(current, tail)?,
    });
  }
  Ok(composed.expect("at least one operand"))
}
