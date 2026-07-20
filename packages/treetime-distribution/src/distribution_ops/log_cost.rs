use crate::policy::{Plain, YAxisPolicy};
use crate::{Distribution, DistributionFunction};
use eyre::Report;
use ndarray::Array1;
use treetime_utils::array::ndarray::min_or;
use treetime_utils::make_error;

/// Multiplies a plain-probability distribution by `exp(-weight(t))`, where the
/// caller supplies `weight` in negative-log space (a cost). The product is
/// formed in negative-log space and shifted by its minimum for numerical
/// stability, then peak-normalized.
///
/// Concrete distributions only: `Formula` is rejected, `Empty` passes through.
pub fn distribution_apply_neg_log_weight<F>(
  distribution: &Distribution<Plain>,
  weight: F,
) -> Result<Distribution<Plain>, Report>
where
  F: Fn(f64) -> Result<f64, Report>,
{
  if matches!(distribution, Distribution::Empty) {
    return Ok(Distribution::Empty);
  }
  if matches!(distribution, Distribution::Formula(_)) {
    return make_error!(
      "distribution_apply_neg_log_weight requires a concrete Point, Range, or Function distribution"
    );
  }

  let times = distribution.t();
  let amplitudes = distribution.y();

  let weights: Array1<f64> = times.iter().map(|&t| weight(t)).collect::<Result<_, Report>>()?;
  // -inf for non-positive amplitudes: makes neg_log +inf, zeroing them after exp.
  let neg_log = weights - amplitudes.mapv(|a| if Plain::is_defined(a) { a.ln() } else { f64::NEG_INFINITY });
  let minimum = min_or(&neg_log, f64::INFINITY);
  if !minimum.is_finite() {
    return make_error!("distribution_apply_neg_log_weight found no finite weight over the distribution grid");
  }
  let scaled = neg_log.mapv(|value| (minimum - value).exp());

  let result = if let Distribution::Function(f) = distribution {
    Distribution::Function(DistributionFunction::from_start_dx_values(f.x_min(), f.dx(), scaled)?)
  } else {
    Distribution::function(times, scaled)?
  };

  Ok(result.normalize())
}
