use crate::Distribution;
use crate::distribution_ops::convolve::distribution_convolution_fine;
use crate::distribution_ops::mass_domain::{
  mass_bounded_domain, peak_normalized_if_mass_sizable, resample_to_mass_window,
};
use crate::policy::NegLog;
use eyre::Report;
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, Side, SoftTailLaw};

/// Convolve two neg-log distributions across a tree edge onto a mass-based output grid.
///
/// The output grid is chosen from the operands before convolving. A convolution's support is the
/// Minkowski sum of the operand supports, and its moments add (mean + mean, variance + variance), so
/// its mass lies within the Minkowski sum of the operand mass-bounded domains: `[lo_a + lo_b,
/// hi_a + hi_b]` bounds the output mass window without the result in hand. That window is conservative
/// (a hair wider than the finished result's own window, since summed operand extents exceed the
/// combined spread) and so never clips mass. The convolution then runs at the finest operand
/// resolution and the result is landed on the window in a single resample, with no intermediate
/// operand-spacing coarsening.
///
/// `soft` names the side that continues past the grid edge with a fitted log-linear tail: a backward
/// message leaves the parent's age unbounded to the past (`Side::Left`), a forward message leaves the
/// node's time unbounded to the future (`Side::Right`). The opposite side is `Hard`, the exact
/// Minkowski support bound where probability is zero beyond; the window's hard edge is clamped to the
/// convolution's own support so it lands exactly there. The soft tail is refit on the landed grid.
///
/// A `Point`, `Range`, or `Empty` convolution result is exact and already compact, so it is returned
/// peak-normalized without windowing. When an operand has no mass-bounded domain (a `Point`/`Range`
/// operand, or a non-mass-sizable tail), the window is sized from the result's own mass instead.
pub fn convolve_across_edge(
  a: &Distribution<NegLog>,
  b: &Distribution<NegLog>,
  soft: Side,
  eps: f64,
  grid_points: usize,
) -> Result<Distribution<NegLog>, Report> {
  let conv = distribution_convolution_fine(a, b)?;
  let Distribution::Function(conv) = conv else {
    return Ok(conv.normalize());
  };
  // keep grid as is (truncate to machine precision range), but refit tails.


  // Fit the soft tail from the fine convolution grid so the window's extrapolated outer points are
  // sensible; declare the opposite side hard. The soft slope is refit once the result is landed.
  let soft_law = SoftTailLaw::fit(conv.grid_fn(), soft, DEFAULT_TAIL_FIT_POINTS)?;
  let conv = match soft {
    Side::Left => conv
      .with_left_extrap(BoundaryBehavior::Linear(soft_law))?
      .with_right_extrap(BoundaryBehavior::Hard)?,
    Side::Right => conv
      .with_right_extrap(BoundaryBehavior::Linear(soft_law))?
      .with_left_extrap(BoundaryBehavior::Hard)?,
  };

  let Some(normalized) = peak_normalized_if_mass_sizable(&conv) else {
    return Ok(Distribution::Function(conv).normalize());
  };

  let (lo, hi) = match convolution_output_window(a, b, eps) {
    Some((mut lo, mut hi)) => {
      // Clamp the hard side to the convolution's own support so the hard edge lands exactly on the
      // zero-beyond bound (never past a cropped FFT tail, which would read `+inf`); the soft side
      // keeps its operand-derived, conservative extent.
      match soft {
        Side::Left => hi = hi.min(normalized.x_max()),
        Side::Right => lo = lo.max(normalized.x_min()),
      }
      if hi > lo {
        (lo, hi)
      } else {
        mass_bounded_domain(&normalized, eps)?
      }
    },
    None => mass_bounded_domain(&normalized, eps)?,
  };

  resample_to_mass_window(&normalized, lo, hi, grid_points)
}

/// The convolution's output mass window from the operands, before convolving.
///
/// Each operand contributes its own mass-bounded domain (a `Point`/`Range` its exact extent, a
/// `Function` its `eps`-mass domain), and the Minkowski sum of the two bounds the result's mass window
/// because convolution moments add. `None` when either operand has no mass-bounded domain (an `Empty`
/// or `Formula` operand, or a non-mass-sizable `Function`), leaving the caller to size by the result.
fn convolution_output_window(a: &Distribution<NegLog>, b: &Distribution<NegLog>, eps: f64) -> Option<(f64, f64)> {
  let (lo_a, hi_a) = operand_mass_domain(a, eps)?;
  let (lo_b, hi_b) = operand_mass_domain(b, eps)?;
  Some((lo_a + lo_b, hi_a + hi_b))
}

/// One operand's mass-bounded domain: a `Point`'s location, a `Range`'s extent, or a `Function`'s
/// `eps`-mass domain. `None` for an `Empty`, `Formula`, or non-mass-sizable operand.
fn operand_mass_domain(dist: &Distribution<NegLog>, eps: f64) -> Option<(f64, f64)> {
  match dist {
    Distribution::Point(p) => Some((p.t(), p.t())),
    Distribution::Range(r) => Some((r.start(), r.end())),
    Distribution::Function(f) => mass_bounded_domain(f, eps).ok(),
    Distribution::Empty | Distribution::Formula(_) => None,
  }
}
