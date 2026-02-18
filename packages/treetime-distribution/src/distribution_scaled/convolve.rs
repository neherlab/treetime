use crate::distribution_ops::convolve::distribution_convolution;
use crate::distribution_scaled::distribution_scaled::ScaledDistribution;
use eyre::Report;

/// Convolve two scaled distributions.
///
/// Delegates to distribution_convolution for the inner (normalized) distributions,
/// then combines log_scales.
pub fn scaled_distribution_convolution(
  a: &ScaledDistribution,
  b: &ScaledDistribution,
) -> Result<ScaledDistribution, Report> {
  if a.is_empty() || b.is_empty() {
    return Ok(ScaledDistribution::default());
  }

  let conv_inner = distribution_convolution(a.inner(), b.inner())?;

  let max_val = conv_inner.max_value();
  if max_val <= 0.0 || !max_val.is_finite() {
    return Ok(ScaledDistribution::default());
  }

  let combined_log_scale = a.log_scale() + b.log_scale() + max_val.ln();
  let normalized_conv = conv_inner.normalize();

  Ok(ScaledDistribution::from_parts(combined_log_scale, normalized_conv))
}
