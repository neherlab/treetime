use crate::policy::{NegLog, Plain, YAxisPolicy};
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
///
/// Takes a closure rather than a pre-evaluated weight array because the grid is
/// not known until variant dispatch. Could take `&Array1<f64>` instead and push
/// dispatch to the caller.
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
    return make_error!("distribution_apply_neg_log_weight requires a concrete Point, Range, or Function distribution");
  }

  let times = distribution.t();
  let amplitudes = distribution.y();

  let weights: Array1<f64> = times.iter().map(|&t| weight(t)).collect::<Result<_, Report>>()?;
  // -inf for non-positive amplitudes: makes neg_log +inf, zeroing them after exp.
  let neg_log = weights
    - amplitudes.mapv(|a| {
      if Plain::is_defined(a) {
        a.ln()
      } else {
        f64::NEG_INFINITY
      }
    });
  let minimum = min_or(&neg_log, f64::INFINITY);
  if !minimum.is_finite() {
    return make_error!("distribution_apply_neg_log_weight found no finite weight over the distribution grid");
  }
  let scaled = neg_log.mapv(|value| (minimum - value).exp());

  let result: Distribution<Plain> = if let Distribution::Function(f) = distribution {
    Distribution::Function(DistributionFunction::from_start_dx_values(f.x_min(), f.dx(), scaled)?)
  } else {
    Distribution::function(times, scaled)?
  };

  Ok(result.normalize())
}

/// Multiplies a negative-log distribution by `exp(-weight(t))`, where the caller supplies
/// `weight` in negative-log space (a cost).
///
/// This is the log-space counterpart of [`distribution_apply_neg_log_weight`]. Under [`NegLog`]
/// the stored ordinate is already `-ln(p)`, so multiplying by `exp(-weight)` is adding the weight
/// to the ordinate: `-ln(p) + cost = -ln(p * exp(-cost))`. No exp round-trip is needed. The result
/// is peak-normalized by subtracting its minimum ordinate.
///
/// Concrete distributions only: `Formula` is rejected, `Empty` passes through.
pub fn distribution_add_neg_log_weight<F>(
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
    return make_error!("distribution_add_neg_log_weight requires a concrete Point, Range, or Function distribution");
  }

  let times = distribution.t();
  let ordinates = distribution.y();

  let weights: Array1<f64> = times.iter().map(|&t| weight(t)).collect::<Result<_, Report>>()?;
  let combined = ordinates + weights;

  let minimum = min_or(&combined, f64::INFINITY);
  if !minimum.is_finite() {
    return make_error!("distribution_add_neg_log_weight found no finite weight over the distribution grid");
  }

  let result: Distribution<NegLog> = if let Distribution::Function(f) = distribution {
    Distribution::Function(DistributionFunction::from_start_dx_values(f.x_min(), f.dx(), combined)?)
  } else {
    Distribution::function(times, combined)?
  };

  Ok(result.normalize())
}
