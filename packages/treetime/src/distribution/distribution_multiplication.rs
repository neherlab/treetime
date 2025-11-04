use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use eyre::Report;
use ndarray::Array1;

/// Multiplies two distributions pointwise (intersection of constraints).
pub fn distribution_multiplication(a: &Distribution, b: &Distribution) -> Result<Distribution, Report> {
  match (a, b) {
    (Distribution::Empty, _) | (_, Distribution::Empty) => Ok(Distribution::empty()),
    (Distribution::Point(a), Distribution::Point(b)) => multiply_point_point(a, b),
    (Distribution::Point(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Point(a)) => {
      multiply_point_function(a, b)
    },
    (Distribution::Point(a), Distribution::Range(b)) | (Distribution::Range(b), Distribution::Point(a)) => {
      multiply_point_range(a, b)
    },
    (Distribution::Range(a), Distribution::Range(b)) => multiply_range_range(a, b),
    (Distribution::Range(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Range(a)) => {
      multiply_range_function(a, b)
    },
    (Distribution::Function(a), Distribution::Function(b)) => multiply_function_function(a, b),
  }
}

fn multiply_point_point(a: &DistributionPoint<f64>, b: &DistributionPoint<f64>) -> Result<Distribution, Report> {
  const EPS: f64 = 1e-9;
  if (a.t() - b.t()).abs() > EPS {
    return Ok(Distribution::empty());
  }
  let amplitude = a.amplitude() * b.amplitude();
  Ok(Distribution::point(a.t(), amplitude))
}

fn multiply_point_range(
  point: &DistributionPoint<f64>,
  range: &crate::distribution::distribution_range::DistributionRange<f64>,
) -> Result<Distribution, Report> {
  const EPS: f64 = 1e-9;
  let t = point.t();
  if t < range.start() - EPS || t > range.end() + EPS {
    return Ok(Distribution::empty());
  }
  let amplitude = point.amplitude() * range.amplitude();
  Ok(Distribution::point(t, amplitude))
}

fn multiply_point_function(
  point: &DistributionPoint<f64>,
  func: &DistributionFunction<f64>,
) -> Result<Distribution, Report> {
  let t = point.t();
  let func_value = func.interp(t).unwrap_or(0.0);
  let amplitude = point.amplitude() * func_value;
  if amplitude <= 0.0 {
    return Ok(Distribution::empty());
  }
  Ok(Distribution::point(t, amplitude))
}

fn multiply_range_range(
  a: &crate::distribution::distribution_range::DistributionRange<f64>,
  b: &crate::distribution::distribution_range::DistributionRange<f64>,
) -> Result<Distribution, Report> {
  let overlap_start = a.start().max(b.start());
  let overlap_end = a.end().min(b.end());

  if overlap_start >= overlap_end {
    return Ok(Distribution::empty());
  }

  let amplitude = a.amplitude() * b.amplitude();
  Ok(Distribution::range((overlap_start, overlap_end), amplitude))
}

fn multiply_range_function(
  range: &crate::distribution::distribution_range::DistributionRange<f64>,
  func: &DistributionFunction<f64>,
) -> Result<Distribution, Report> {
  const EPS: f64 = 1e-9;

  let func_t = func.t();
  let func_y = func.y();

  // Filter to range support
  let filtered: Vec<(f64, f64)> = func_t
    .iter()
    .zip(func_y.iter())
    .filter(|&(&t, _)| t >= range.start() - EPS && t <= range.end() + EPS)
    .map(|(&t, &y)| (t, y * range.amplitude()))
    .collect();

  if filtered.is_empty() {
    return Ok(Distribution::empty());
  }

  let times: Vec<f64> = filtered.iter().map(|(t, _)| *t).collect();
  let values: Vec<f64> = filtered.iter().map(|(_, v)| *v).collect();
  Distribution::function(Array1::from_vec(times), Array1::from_vec(values))
}

fn multiply_function_function(
  a: &DistributionFunction<f64>,
  b: &DistributionFunction<f64>,
) -> Result<Distribution, Report> {
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

  // Create uniform grid in overlap region
  let n_points = a.t().len().max(b.t().len());
  let times = Array1::linspace(overlap_min, overlap_max, n_points);

  let values: Array1<f64> = times
    .iter()
    .map(|&t| {
      let va = a.interp(t).unwrap_or(0.0);
      let vb = b.interp(t).unwrap_or(0.0);
      va * vb
    })
    .collect();

  Distribution::function(times, values)
}
