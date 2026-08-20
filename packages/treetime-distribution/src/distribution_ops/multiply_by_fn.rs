use crate::distribution_ops::mass_domain::refit_soft_tails;
use crate::policy::NegLog;
use crate::{Distribution, DistributionFunction};
use eyre::Report;
use ndarray::Array1;
use treetime_utils::array::ndarray::min_or;
use treetime_utils::make_error;

/// Multiplies a negative-log distribution by `exp(-weight(t))`, where the caller supplies
/// `weight` in negative-log space (a cost).
///
/// Under [`NegLog`] the stored ordinate is already `-ln(p)`, so multiplying by `exp(-weight)` is
/// adding the weight to the ordinate: `-ln(p) + cost = -ln(p * exp(-cost))`. No exp round-trip is
/// needed. The result is peak-normalized by subtracting its minimum ordinate.
///
/// Concrete distributions only: `Formula` is rejected, `Empty` passes through.
///
/// Tail policy survives the addition (`kb/decisions/distribution-tails-and-arithmetic.md`). Adding a
/// finite per-point cost is a pointwise ordinate operation: it can neither turn zero probability
/// positive nor make an undefined region defined, so each side keeps the input distribution's tail
/// *class* -- a hard or hard-approach side stays hard, an `Error` side stays `Error`, a soft side
/// stays soft. Unlike `normalize`'s constant shift the cost varies with `t`, so it changes the soft
/// tail's local slope; the `Linear` law is therefore re-fit from the combined grid, while the
/// edge-relative `Hard`/`HardApproach` laws carry through unchanged. Erasing the policy to `Error`
/// here disabled the mass re-window (design D3) for every coalescent internal and root node, since
/// that node's contribution routes through this function in the backward pass.
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
  let ordinates = distribution.y();

  let weights: Array1<f64> = times.iter().map(|&t| weight(t)).collect::<Result<_, Report>>()?;
  let combined = ordinates + weights;

  let minimum = min_or(&combined, f64::INFINITY);
  if !minimum.is_finite() {
    return make_error!("distribution_multiply_by_fn found no finite weight over the distribution grid");
  }

  let result: Distribution<NegLog> = if let Distribution::Function(f) = distribution {
    // Rebuild carries no tails, so re-declare the input's per-side policy, then re-fit the soft slope
    // the varying weight changed (Hard/HardApproach/Error are carried unchanged).
    let combined_fn = DistributionFunction::from_start_dx_values(f.x_min(), f.dx(), combined)?
      .with_left_extrap(f.left_extrap())?
      .with_right_extrap(f.right_extrap())?;
    Distribution::Function(refit_soft_tails(combined_fn)?)
  } else {
    // Point/Range carry no interpolated tail (`Error` on both sides), so the fresh function is correct.
    Distribution::function(times, combined)?
  };

  Ok(result.normalize())
}
