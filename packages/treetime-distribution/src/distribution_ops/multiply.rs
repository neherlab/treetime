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
use ndarray::Array1;

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
  let func_value = func.interp(t);
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
  match distribution_support_intersection((range.start(), range.end()), (func.x_min(), func.x_max())) {
    SupportIntersection::Disjoint => Ok(Distribution::empty()),
    SupportIntersection::Point(t) => {
      let amplitude = Y::multiply(range.amplitude(), func.interp(t));
      Ok(Distribution::point(t, amplitude))
    },
    SupportIntersection::Interval(bounds) => {
      let n_points = distribution_support_n_points(bounds, func.dx())?;
      let function =
        DistributionFunction::from_n_points(bounds, n_points, |t| Y::multiply(range.amplitude(), func.interp(t)))?;
      Ok(Distribution::Function(function))
    },
  }
}

fn multiply_function_function<Y: YAxisPolicy>(
  a: &DistributionFunction<f64, Y>,
  b: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  match distribution_support_intersection((a.x_min(), a.x_max()), (b.x_min(), b.x_max())) {
    SupportIntersection::Disjoint => Ok(Distribution::empty()),
    SupportIntersection::Point(t) => Ok(Distribution::point(t, Y::multiply(a.interp(t), b.interp(t)))),
    SupportIntersection::Interval(bounds) => {
      let function =
        DistributionFunction::from_n_points(bounds, a.len().max(b.len()), |t| Y::multiply(a.interp(t), b.interp(t)))?;
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
  let overlap_min = a.t_min().max(b.x_min());
  let overlap_max = a.t_max().min(b.x_max());

  if overlap_min >= overlap_max {
    return Ok(Distribution::empty());
  }

  // Evaluate the Formula on the Function's grid (matching multiply_function_function).
  // The Function provides the grid; the Formula is evaluated at each grid point.
  let n_points = b.len();
  let values = (0..n_points)
    .map(|i| {
      let t = overlap_min + (overlap_max - overlap_min) * (i as f64 / (n_points - 1) as f64);
      let va = a.eval_single(t)?;
      let vb = b.interp(t);
      Ok(Y::multiply(va, vb))
    })
    .collect::<Result<Vec<f64>, Report>>()?;

  let distribution_fn = DistributionFunction::from_range_values((overlap_min, overlap_max), Array1::from_vec(values))?;
  Ok(Distribution::Function(distribution_fn))
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
