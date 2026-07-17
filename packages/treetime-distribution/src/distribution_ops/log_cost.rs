use crate::Distribution;
use crate::policy::Plain;
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
    Distribution::Point(_) | Distribution::Range(_) | Distribution::Function(_) => {
      let times = distribution.t();
      let amplitudes = distribution.y();
      let neg_log = amplitudes
        .iter()
        .zip(times.iter())
        .map(|(&amplitude, &time)| Ok::<f64, Report>(weight(time)? - amplitude.ln()))
        .collect::<Result<Vec<_>, Report>>()?;
      let Some(minimum) = neg_log
        .iter()
        .copied()
        .reduce(f64::min)
        .filter(|minimum| minimum.is_finite())
      else {
        return make_error!("distribution_apply_neg_log_weight found no finite weight over the distribution grid");
      };
      let scaled = Array1::from_iter(neg_log.iter().map(|&value| (minimum - value).exp()));
      Ok(Distribution::function(times, scaled)?.normalize())
    },
  }
}
