use crate::Distribution;
use crate::distribution_formula::DistributionFormula;
use crate::distribution_function::DistributionFunction;
use crate::distribution_point::DistributionPoint;
use crate::distribution_range::DistributionRange;
use crate::y_axis_policy::YAxisPolicy;
use eyre::Report;
use ndarray::Array1;
use std::sync::Arc;

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
  const EPS: f64 = 1e-9;

  let func_t = func.t();
  let func_y = func.y();

  // Filter to range support
  let filtered: Vec<(f64, f64)> = func_t
    .iter()
    .zip(func_y.iter())
    .filter(|&(&t, _)| t >= range.start() - EPS && t <= range.end() + EPS)
    .map(|(&t, &y)| (t, Y::multiply(y, range.amplitude())))
    .collect();

  if filtered.is_empty() {
    return Ok(Distribution::empty());
  }

  let overlap_min = filtered.first().unwrap().0;
  let overlap_max = filtered.last().unwrap().0;
  let values: Vec<f64> = filtered.iter().map(|(_, v)| *v).collect();
  let values_array = Array1::from_vec(values);

  let distribution_fn = DistributionFunction::from_range_values((overlap_min, overlap_max), values_array)?;
  Ok(Distribution::Function(distribution_fn))
}

fn multiply_function_function<Y: YAxisPolicy>(
  a: &DistributionFunction<f64, Y>,
  b: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  // Find overlapping support
  let a_min = a.t().first().copied().unwrap_or(0.0);
  let a_max = a.t().last().copied().unwrap_or(0.0);
  let b_min = b.t().first().copied().unwrap_or(0.0);
  let b_max = b.t().last().copied().unwrap_or(0.0);

  let overlap_min = a_min.max(b_min);
  let overlap_max = a_max.min(b_max);

  if overlap_min >= overlap_max {
    return Ok(Distribution::empty());
  }

  let n_points = a.t().len().max(b.t().len());

  let values: Array1<f64> = Array1::from_shape_fn(n_points, |i| {
    let t = overlap_min + (overlap_max - overlap_min) * (i as f64 / (n_points - 1) as f64);
    let va = a.interp(t);
    let vb = b.interp(t);
    Y::multiply(va, vb)
  });

  let distribution_fn = DistributionFunction::from_range_values((overlap_min, overlap_max), values)?;
  Ok(Distribution::Function(distribution_fn))
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

  let a = Arc::new(a.clone());
  let b = Arc::new(b.clone());

  let eval_fn = move |t: f64| -> eyre::Result<f64> {
    let va = a.eval_single(t)?;
    let vb = b.eval_single(t)?;
    Ok(Y::multiply(va, vb))
  };

  Ok(Distribution::Formula(DistributionFormula::new(
    eval_fn,
    overlap_min,
    overlap_max,
  )))
}

fn multiply_formula_function<Y: YAxisPolicy>(
  a: &DistributionFormula<Y>,
  b: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let b_min = b.t().first().copied().unwrap_or(0.0);
  let b_max = b.t().last().copied().unwrap_or(0.0);

  let overlap_min = a.t_min().max(b_min);
  let overlap_max = a.t_max().min(b_max);

  if overlap_min >= overlap_max {
    return Ok(Distribution::empty());
  }

  let a = Arc::new(a.clone());
  let b = Arc::new(b.clone());

  let eval_fn = move |t: f64| -> eyre::Result<f64> {
    let va = a.eval_single(t)?;
    let vb = b.interp(t);
    Ok(Y::multiply(va, vb))
  };

  Ok(Distribution::Formula(DistributionFormula::new(
    eval_fn,
    overlap_min,
    overlap_max,
  )))
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

  let a = Arc::new(a.clone());
  let b_amplitude = b.amplitude();

  let eval_fn = move |t: f64| -> eyre::Result<f64> {
    let va = a.eval_single(t)?;
    Ok(Y::multiply(va, b_amplitude))
  };

  Ok(Distribution::Formula(DistributionFormula::new(
    eval_fn,
    overlap_min,
    overlap_max,
  )))
}
