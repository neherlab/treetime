use crate::Distribution;
use crate::policy::Plain;
use eyre::Report;
use ndarray::Array1;
use treetime_utils::make_error;

/// Multiplies a plain-probability distribution by `exp(-weight(t))`, where the
/// caller supplies `weight` in negative-log space (a cost). The weight is
/// evaluated on the distribution's own grid and the result is peak-normalized.
///
/// The product is formed with a log-sum-exp shift for numerical stability: each
/// amplitude is scaled by `exp(min_j weight_j - weight_i)`, so the smallest
/// weight maps to a unit factor and the exponential cannot overflow. The final
/// `normalize()` divides by the peak.
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
      let weights = times
        .iter()
        .map(|&time| weight(time))
        .collect::<Result<Vec<_>, Report>>()?;
      let Some(minimum) = weights.iter().copied().reduce(f64::min).filter(|minimum| minimum.is_finite()) else {
        return make_error!("distribution_apply_neg_log_weight found no finite weight over the distribution grid");
      };
      let amplitudes = Array1::from_iter(
        distribution
          .y()
          .iter()
          .zip(weights)
          .map(|(&amplitude, weight)| amplitude * (minimum - weight).exp()),
      );
      Ok(Distribution::function(times, amplitudes)?.normalize())
    },
  }
}
