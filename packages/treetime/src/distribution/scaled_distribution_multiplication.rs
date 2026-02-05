use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_multiplication::distribution_multiplication;
use crate::distribution::scaled_distribution::ScaledDistribution;
use crate::distribution::y_axis_policy::Plain;
use approx::ulps_eq;
use eyre::Report;
use itertools::Itertools;
use ndarray::Array1;
use treetime_ops::multiply_many_lazy_normalize;

struct AlignedFunctionArrays<'a> {
  arrays: Vec<&'a Array1<f64>>,
  x_min: f64,
  dx: f64,
}

fn try_extract_aligned_function_arrays<'a>(
  distributions: &'a [&'a ScaledDistribution],
) -> Option<AlignedFunctionArrays<'a>> {
  let first = distributions.first()?;
  let Distribution::Function(first_func) = first.inner() else {
    return None;
  };

  let x_min = first_func.x_min();
  let dx = first_func.dx();
  let len = first_func.len();

  let mut arrays = Vec::with_capacity(distributions.len());

  for dist in distributions {
    match dist.inner() {
      Distribution::Function(f) => {
        if !ulps_eq!(f.x_min(), x_min, max_ulps = 10) || !ulps_eq!(f.dx(), dx, max_ulps = 10) || f.len() != len {
          return None;
        }
        arrays.push(f.y());
      },
      _ => return None,
    }
  }

  Some(AlignedFunctionArrays { arrays, x_min, dx })
}

/// Multiply two scaled distributions.
///
/// Delegates to distribution_multiplication for the inner (normalized) distributions,
/// then combines log_scales.
pub fn scaled_distribution_multiplication(
  a: &ScaledDistribution,
  b: &ScaledDistribution,
) -> Result<ScaledDistribution, Report> {
  if a.is_empty() || b.is_empty() {
    return Ok(ScaledDistribution::default());
  }

  let product_inner = distribution_multiplication(a.inner(), b.inner())?;

  let max_val = product_inner.max_value();
  if max_val <= 0.0 || !max_val.is_finite() {
    return Ok(ScaledDistribution::default());
  }

  let combined_log_scale = a.log_scale() + b.log_scale() + max_val.ln();
  let normalized_product = product_inner.normalize();

  Ok(ScaledDistribution::from_parts(combined_log_scale, normalized_product))
}

/// Multiply many scaled distributions.
///
/// When all inputs are Function distributions with matching grids, delegates
/// to `treetime_ops::multiply_many_lazy_normalize` for validated numerical
/// stability. Otherwise falls back to pairwise multiplication.
pub fn scaled_distribution_multiply_many(distributions: &[&ScaledDistribution]) -> Result<ScaledDistribution, Report> {
  match distributions {
    [] => Ok(ScaledDistribution::default()),
    [single] => Ok((*single).clone()),
    _ => {
      if distributions.iter().any(|d| d.is_empty()) {
        return Ok(ScaledDistribution::default());
      }

      let total_input_log_scale: f64 = distributions.iter().map(|d| d.log_scale()).sum();

      if let Some(aligned) = try_extract_aligned_function_arrays(distributions) {
        let array_refs = aligned.arrays.iter().copied().collect_vec();
        let result = multiply_many_lazy_normalize(&array_refs);

        if !result.log_scale.is_finite() {
          return Ok(ScaledDistribution::default());
        }

        let result_func =
          DistributionFunction::<f64, Plain>::from_start_dx_values(aligned.x_min, aligned.dx, result.normalized)?;
        let total_log_scale = total_input_log_scale + result.log_scale;

        return Ok(ScaledDistribution::from_parts(
          total_log_scale,
          Distribution::Function(result_func),
        ));
      }

      let mut product = distributions[0].inner().clone();
      for dist in distributions.iter().skip(1) {
        product = distribution_multiplication(&product, dist.inner())?;
      }

      let max_val = product.max_value();
      if max_val <= 0.0 || !max_val.is_finite() {
        return Ok(ScaledDistribution::default());
      }

      Ok(ScaledDistribution::from_parts(
        total_input_log_scale + max_val.ln(),
        product.normalize(),
      ))
    },
  }
}
