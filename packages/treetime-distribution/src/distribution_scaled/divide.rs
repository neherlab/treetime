use crate::distribution_ops::divide::distribution_division;
use crate::distribution_scaled::scaled::ScaledDistribution;
use eyre::Report;
use treetime_utils::make_error;

/// Divide one scaled distribution by another.
///
/// Delegates to distribution_division for the inner (normalized) distributions,
/// then combines log_scales (subtraction for division).
pub fn scaled_distribution_division(
  dividend: &ScaledDistribution,
  divisor: &ScaledDistribution,
) -> Result<ScaledDistribution, Report> {
  if dividend.is_empty() {
    return Ok(ScaledDistribution::default());
  }

  if divisor.is_empty() {
    return make_error!("Cannot divide by empty distribution");
  }

  let quotient_inner = distribution_division(dividend.inner(), divisor.inner())?;

  let max_val = quotient_inner.max_value();
  if max_val <= 0.0 || !max_val.is_finite() {
    return Ok(ScaledDistribution::default());
  }

  let combined_log_scale = dividend.log_scale() - divisor.log_scale() + max_val.ln();
  let normalized_quotient = quotient_inner.normalize();

  Ok(ScaledDistribution::from_parts(combined_log_scale, normalized_quotient))
}
