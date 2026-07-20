use crate::policy::{Plain, YAxisPolicy};
use crate::{Distribution, DistributionFunction};
use eyre::Report;
use ndarray::Array1;
use treetime_utils::make_error;

/// Multiplies a plain-probability distribution by `exp(-weight(t))`, where the
/// caller supplies `weight` in negative-log space (a cost). The weight is
/// evaluated on the distribution's own grid and the result is peak-normalized.
///
/// The product is formed in negative-log space for numerical stability. Each
/// grid point's combined value is `weight_i - ln(amplitude_i)`; amplitudes are
/// then `exp(min_j(combined_j) - combined_i)`. Shifting by the minimum of the
/// *combined* value, rather than the minimum of the weight alone, places the
/// product peak at unit magnitude. A peak that sits at a high-weight grid point
/// therefore survives, instead of underflowing `exp(min_weight - weight_i)` to
/// zero when the weight spans more than ~700 across the grid.
///
/// A grid point whose amplitude is not a valid plain probability (`amplitude <=
/// 0`, per `Plain::is_defined`) carries no mass: its combined value is treated
/// as `+inf`, excluding it from the shift minimum and mapping it to `0` after
/// exponentiation. This matches how `distribution_multiplication` already drops
/// non-positive amplitudes, and guards `ln(amplitude)` against a slightly
/// negative amplitude produced by convolution/interpolation undershoot near a
/// grid tail (order `-1e-26`), which would otherwise be `NaN` and collapse the
/// whole distribution to `Empty` when normalized.
///
/// Concrete distributions only: a `Formula` has no grid to evaluate on and is
/// rejected. `Empty` passes through unchanged.
pub fn distribution_apply_neg_log_weight<F>(
  distribution: &Distribution<Plain>,
  weight: F,
) -> Result<Distribution<Plain>, Report>
where
  F: Fn(f64) -> Result<f64, Report>,
{
  match distribution {
    Distribution::Empty => Ok(Distribution::Empty),
    Distribution::Formula(_) => {
      make_error!("distribution_apply_neg_log_weight requires a concrete Point, Range, or Function distribution")
    },
    Distribution::Point(_) | Distribution::Range(_) => {
      let times = distribution.t();
      let amplitudes = distribution.y();
      let scaled = apply_neg_log_weight(&times, &amplitudes, &weight)?;
      Ok(Distribution::function(times, scaled)?.normalize())
    },
    Distribution::Function(function) => {
      let times = function.t();
      let scaled = apply_neg_log_weight(&times, function.y(), &weight)?;
      let function = DistributionFunction::from_start_dx_values(function.x_min(), function.dx(), scaled)?;
      Ok(Distribution::Function(function).normalize())
    },
  }
}

fn apply_neg_log_weight<F>(times: &Array1<f64>, amplitudes: &Array1<f64>, weight: &F) -> Result<Array1<f64>, Report>
where
  F: Fn(f64) -> Result<f64, Report>,
{
  let neg_log = amplitudes
    .iter()
    .zip(times.iter())
    .map(|(&amplitude, &time)| {
      // Non-positive amplitude = no mass. Map to +inf cost instead of taking
      // `ln` (negative -> NaN, zero -> -inf), so the point drops out cleanly.
      if Plain::is_defined(amplitude) {
        Ok::<f64, Report>(weight(time)? - amplitude.ln())
      } else {
        Ok(f64::INFINITY)
      }
    })
    .collect::<Result<Vec<_>, Report>>()?;
  let Some(minimum) = neg_log
    .iter()
    .copied()
    .reduce(f64::min)
    .filter(|minimum| minimum.is_finite())
  else {
    return make_error!("distribution_apply_neg_log_weight found no finite weight over the distribution grid");
  };
  Ok(Array1::from_iter(neg_log.iter().map(|&value| (minimum - value).exp())))
}
