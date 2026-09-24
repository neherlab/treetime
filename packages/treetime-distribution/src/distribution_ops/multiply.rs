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

const FORMULA_GRID_SIZE: usize = 200;

pub fn distribution_multiplication<Y: YAxisPolicy>(
  a: &Distribution<Y>,
  b: &Distribution<Y>,
) -> Result<Distribution<Y>, Report> {
  match (a, b) {
    (Distribution::Empty, _) | (_, Distribution::Empty) => {
      let a_domain = distribution_hard_domain(a);
      let b_domain = distribution_hard_domain(b);
      guarded_empty_result("multiplication", a_domain, b_domain)
    },
    (Distribution::Point(a), Distribution::Point(b)) => multiply_point_point::<Y>(a, b),
    (Distribution::Point(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Point(a)) => {
      multiply_point_function::<Y>(a, b)
    },
    (Distribution::Point(a), Distribution::Range(b)) | (Distribution::Range(b), Distribution::Point(a)) => {
      multiply_point_range::<Y>(a, b)
    },
    (Distribution::Range(a), Distribution::Range(b)) => multiply_range_range::<Y>(a, b),
    (Distribution::Range(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Range(a)) => {
      multiply_range_function::<Y>(a, b)
    },
    (Distribution::Function(a), Distribution::Function(b)) => multiply_function_function::<Y>(a, b),
    (Distribution::Formula(a), Distribution::Formula(b)) => multiply_formula_formula::<Y>(a, b),
    (Distribution::Formula(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Formula(a)) => {
      multiply_formula_function::<Y>(a, b)
    },
    (Distribution::Formula(a), Distribution::Point(b)) | (Distribution::Point(b), Distribution::Formula(a)) => {
      multiply_formula_point::<Y>(a, b)
    },
    (Distribution::Formula(a), Distribution::Range(b)) | (Distribution::Range(b), Distribution::Formula(a)) => {
      multiply_formula_range::<Y>(a, b)
    },
  }
}

fn multiply_point_point<Y: YAxisPolicy>(
  a: &DistributionPoint<f64, Y>,
  b: &DistributionPoint<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  const EPS: f64 = 1e-9;
  if (a.t() - b.t()).abs() > EPS {
    let a_domain = Some(point_hard_domain(a));
    let b_domain = Some(point_hard_domain(b));
    return guarded_empty_result("multiplication", a_domain, b_domain);
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
    let point_domain = Some(point_hard_domain(point));
    let range_domain = Some(range_hard_domain(range));
    return guarded_empty_result("multiplication", point_domain, range_domain);
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
    if !tail.is_soft() {
      let point_domain = Some(point_hard_domain(point));
      let function_domain = Some(function_hard_domain(func));
      return guarded_empty_result("multiplication", point_domain, function_domain);
    }
  }
  let func_value = func.interp(t)?;
  let amplitude = Y::multiply(point.amplitude(), func_value);
  if !Y::is_defined(amplitude) {
    let point_domain = Some(point_hard_domain(point));
    let function_domain = Some(function_hard_domain(func));
    return guarded_empty_result("multiplication", point_domain, function_domain);
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
    let a_domain = Some(range_hard_domain(a));
    let b_domain = Some(range_hard_domain(b));
    return guarded_empty_result("multiplication", a_domain, b_domain);
  }

  let amplitude = Y::multiply(a.amplitude(), b.amplitude());
  Ok(Distribution::range((overlap_start, overlap_end), amplitude))
}

fn multiply_range_function<Y: YAxisPolicy>(
  range: &DistributionRange<f64, Y>,
  func: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let a_bounds = (range.start(), range.end());
  let a_tails = (BoundaryBehavior::Hard, BoundaryBehavior::Hard);
  let b_bounds = (func.x_min(), func.x_max());
  let b_tails = (func.left_extrap(), func.right_extrap());
  match multiplication_support_intersection(&[(a_bounds, a_tails), (b_bounds, b_tails)]) {
    SupportIntersection::Disjoint => {
      let a_domain = Some((a_bounds, a_tails));
      let b_domain = Some((b_bounds, b_tails));
      guarded_empty_result("multiplication", a_domain, b_domain)
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
  multiply_functions(&[a, b])
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
pub(crate) fn multiply_functions<Y: YAxisPolicy>(
  functions: &[&DistributionFunction<f64, Y>],
) -> Result<Distribution<Y>, Report> {
  let ordered = canonical_operand_order(functions);
  let (&first, rest) = ordered
    .split_first()
    .expect("multiply_functions requires at least one operand");

  let domains: Vec<HardDomain> = ordered.iter().copied().map(function_hard_domain).collect();
  match multiplication_support_intersection(&domains) {
    SupportIntersection::Disjoint => Ok(Distribution::empty()),
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

fn canonical_operand_order<'a, Y: YAxisPolicy>(
  functions: &[&'a DistributionFunction<f64, Y>],
) -> Vec<&'a DistributionFunction<f64, Y>> {
  let mut ordered = functions.to_vec();
  ordered.sort_by_key(|f| (OrderedFloat(f.x_min()), OrderedFloat(f.x_max()), OrderedFloat(f.dx())));
  ordered
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
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
  let a_domain = formula_hard_domain(a);
  let b_domain = formula_hard_domain(b);
  match multiplication_support_intersection(&[a_domain, b_domain]) {
    SupportIntersection::Disjoint => {
      let a_domain = Some(a_domain);
      let b_domain = Some(b_domain);
      guarded_empty_result("multiplication", a_domain, b_domain)
    },
    SupportIntersection::Point(t) => Ok(Distribution::point(
      t,
      Y::multiply(a.eval_single(t)?, b.eval_single(t)?),
    )),
    SupportIntersection::Interval(bounds) => {
      let grid = Array1::linspace(bounds.0, bounds.1, FORMULA_GRID_SIZE);
      let a_values = a.eval_many(&grid)?;
      let b_values = b.eval_many(&grid)?;
      let values = Zip::from(&a_values)
        .and(&b_values)
        .map_collect(|&a_value, &b_value| Y::multiply(a_value, b_value));
      let function = with_composed_tails(
        DistributionFunction::from_range_values(bounds, values)?,
        a_domain.1,
        b_domain.1,
      )?;
      Ok(Distribution::Function(function))
    },
  }
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
      let a_domain = Some((a_bounds, a_tails));
      let b_domain = Some((b_bounds, b_tails));
      guarded_empty_result("multiplication", a_domain, b_domain)
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
  let a_domain = formula_hard_domain(a);
  let b_domain = range_hard_domain(b);
  match multiplication_support_intersection(&[a_domain, b_domain]) {
    SupportIntersection::Disjoint => {
      let a_domain = Some(a_domain);
      let b_domain = Some(b_domain);
      guarded_empty_result("multiplication", a_domain, b_domain)
    },
    SupportIntersection::Point(t) => Ok(Distribution::point(t, Y::multiply(a.eval_single(t)?, b.amplitude()))),
    SupportIntersection::Interval(bounds) => {
      let grid = Array1::linspace(bounds.0, bounds.1, FORMULA_GRID_SIZE);
      let values = a.eval_many(&grid)?.mapv(|value| Y::multiply(value, b.amplitude()));
      let function = with_composed_tails(
        DistributionFunction::from_range_values(bounds, values)?,
        a_domain.1,
        b_domain.1,
      )?;
      Ok(Distribution::Function(function))
    },
  }
}

#[expect(clippy::float_cmp, reason = "equal bounds define a point support exactly")]
pub(super) fn multiplication_support_intersection(domains: &[HardDomain]) -> SupportIntersection {
  let (hard_lo, soft_lo) = side_bounds(domains, Side::Left);
  let (hard_hi, soft_hi) = side_bounds(domains, Side::Right);

  if izip!(hard_lo, hard_hi).any(|(lo, hi)| lo > hi) {
    return SupportIntersection::Disjoint;
  }

  let lo = hard_lo.or(soft_lo).unwrap_or(f64::INFINITY);
  let hi = hard_hi.or(soft_hi).unwrap_or(f64::NEG_INFINITY);
  if lo == hi {
    SupportIntersection::Point(lo)
  } else if lo < hi {
    SupportIntersection::Interval((lo, hi))
  } else {
    SupportIntersection::Disjoint
  }
}

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

fn with_composed_tails<Y: YAxisPolicy>(
  function: DistributionFunction<f64, Y>,
  a_tails: (BoundaryBehavior, BoundaryBehavior),
  b_tails: (BoundaryBehavior, BoundaryBehavior),
) -> Result<DistributionFunction<f64, Y>, Report> {
  function
    .with_left_extrap(compose_multiplication_tail(a_tails.0, b_tails.0)?)?
    .with_right_extrap(compose_multiplication_tail(a_tails.1, b_tails.1)?)
}

fn compose_multiplication_tail(a: BoundaryBehavior, b: BoundaryBehavior) -> Result<BoundaryBehavior, Report> {
  match (a, b) {
    (BoundaryBehavior::Error, _) | (_, BoundaryBehavior::Error) => Ok(BoundaryBehavior::Error),
    (BoundaryBehavior::Hard, BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_))
    | (BoundaryBehavior::HardApproach(_), BoundaryBehavior::Hard) => Ok(BoundaryBehavior::Hard),
    (BoundaryBehavior::HardApproach(_), BoundaryBehavior::HardApproach(_)) => make_internal_error!(
      "Cannot multiply two HardApproach tails: their product is not representable by a \
       single-parameter hard-approach law, and this composition is unreachable in the inference pipeline"
    ),
    (hard @ (BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_)), _) => Ok(hard),
    (_, hard @ (BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_))) => Ok(hard),
    (BoundaryBehavior::Linear(a_law), BoundaryBehavior::Linear(b_law)) => {
      Ok(BoundaryBehavior::Linear(a_law.compose_multiply(&b_law)))
    },
  }
}

pub(crate) fn distribution_hard_domain<Y: YAxisPolicy>(d: &Distribution<Y>) -> Option<HardDomain> {
  match d {
    Distribution::Empty => None,
    Distribution::Point(p) => Some(point_hard_domain(p)),
    Distribution::Range(r) => Some(range_hard_domain(r)),
    Distribution::Formula(f) => Some(formula_hard_domain(f)),
    Distribution::Function(f) => (!f.is_empty()).then(|| function_hard_domain(f)),
  }
}

pub(crate) fn point_hard_domain<Y: YAxisPolicy>(p: &DistributionPoint<f64, Y>) -> HardDomain {
  ((p.t(), p.t()), (BoundaryBehavior::Hard, BoundaryBehavior::Hard))
}

pub(crate) fn range_hard_domain<Y: YAxisPolicy>(r: &DistributionRange<f64, Y>) -> HardDomain {
  ((r.start(), r.end()), (BoundaryBehavior::Hard, BoundaryBehavior::Hard))
}

pub(crate) fn function_hard_domain<Y: YAxisPolicy>(f: &DistributionFunction<f64, Y>) -> HardDomain {
  ((f.x_min(), f.x_max()), (f.left_extrap(), f.right_extrap()))
}

fn formula_hard_domain<Y: YAxisPolicy>(f: &DistributionFormula<Y>) -> HardDomain {
  (
    (f.t_min(), f.t_max()),
    (BoundaryBehavior::Error, BoundaryBehavior::Error),
  )
}

pub(crate) fn guarded_empty_result<Y: YAxisPolicy>(
  operation: &str,
  a: Option<HardDomain>,
  b: Option<HardDomain>,
) -> Result<Distribution<Y>, Report> {
  let disjoint = match (a, b) {
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

pub(crate) type HardDomain = ((f64, f64), (BoundaryBehavior, BoundaryBehavior));

pub(crate) fn hard_domains_disjoint(
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
