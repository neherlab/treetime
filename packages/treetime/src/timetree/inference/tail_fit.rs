use eyre::Report;
use treetime_distribution::{BoundaryBehavior, Distribution, NegLog};
use treetime_grid::{DEFAULT_TAIL_FIT_POINTS, Side, SoftTailLaw};
use treetime_utils::make_internal_report;

/// Fit the log-linear soft tail on one side of a convolved timetree message.
///
/// Each timetree message has one *soft* side that continues past the grid edge: a backward message
/// leaves the parent's age unbounded to the past, a forward message leaves the node's time unbounded
/// to the future. That side must carry an *integrable* tail. A flat [`BoundaryBehavior::Constant`]
/// has infinite mass and corrupts the quantile and HPD integrals; the log-linear [`SoftTailLaw`]
/// decays and has finite mass (kb/decisions/distribution-tails-and-arithmetic.md).
///
/// The slope is fit from the outermost [`DEFAULT_TAIL_FIT_POINTS`] grid points. The FFT convolution
/// has already reconstructed those points by log-linear extrapolation below its roundoff floor, so
/// the fit recovers the message's own decay rather than reading the untrusted FFT tail. A grid too
/// degenerate to fit is an error, never a silent flat fallback: the total `BoundaryBehavior` enum
/// has no law-less `Linear`, and a fabricated flat tail would reintroduce the non-integrable
/// `Constant` behavior this replaces.
///
/// A `Point`, `Range`, or `Empty` message carries no grid tail; the returned `Constant` is inert
/// because [`Distribution::with_left_extrap`] and [`Distribution::with_right_extrap`] are no-ops on
/// non-`Function` variants.
pub fn fit_message_soft_tail(message: &Distribution<NegLog>, side: Side) -> Result<BoundaryBehavior, Report> {
  let Distribution::Function(function) = message else {
    return Ok(BoundaryBehavior::Constant);
  };
  let law = SoftTailLaw::fit(function.grid_fn(), side, DEFAULT_TAIL_FIT_POINTS).ok_or_else(|| {
    make_internal_report!("Timetree message cannot fit a soft tail on the {side:?} side: the convolved grid is degenerate")
  })?;
  Ok(BoundaryBehavior::Linear(law))
}
