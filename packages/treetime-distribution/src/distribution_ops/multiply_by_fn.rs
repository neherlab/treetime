use crate::distribution_ops::mass_domain::refit_soft_tails;
use crate::policy::NegLog;
use crate::{Distribution, DistributionFunction};
use eyre::Report;
use ndarray::Array1;
use treetime_utils::array::ndarray::min_or;
use treetime_utils::make_error;

pub fn distribution_multiply_by_fn<F>(
  distribution: &Distribution<NegLog>,
  weight: F,
) -> Result<Distribution<NegLog>, Report>
where
  F: Fn(f64) -> Result<f64, Report>,
{
  if matches!(distribution, Distribution::Empty) {
    return Ok(Distribution::Empty);
  }
  if matches!(distribution, Distribution::Formula(_)) {
    return make_error!("distribution_multiply_by_fn requires a concrete Point, Range, or Function distribution");
  }

  let times = distribution.t();
  let ordinates = distribution.y()?;

  let weights: Array1<f64> = times.iter().map(|&t| weight(t)).collect::<Result<_, Report>>()?;
  let combined = ordinates + weights;

  let minimum = min_or(&combined, f64::INFINITY);
  if !minimum.is_finite() {
    return make_error!("distribution_multiply_by_fn found no finite weight over the distribution grid");
  }

  let result: Distribution<NegLog> = if let Distribution::Function(f) = distribution {
    let combined_fn = DistributionFunction::from_start_dx_values(f.x_min(), f.dx(), combined)?
      .with_left_extrap(f.left_extrap())?
      .with_right_extrap(f.right_extrap())?;
    Distribution::Function(refit_soft_tails(combined_fn)?)
  } else {
    Distribution::function(times, combined)?
  };

  Ok(result.normalize())
}
