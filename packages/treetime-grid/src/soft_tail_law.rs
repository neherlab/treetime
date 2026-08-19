use crate::hard_approach_law::Side;
use crate::{GridEdge, GridFn};
use eyre::Report;
use serde::{Deserialize, Serialize};
use treetime_utils::least_squares::LineFit;
use treetime_utils::make_error;

/// Log-linear approach for a *soft* grid boundary: an exponential probability tail that
/// continues the distribution past the grid edge.
///
/// A soft boundary is a representation choice, not a fact about the distribution: the grid edge
/// is only where interpolation stopped, and the density continues beyond it. In negative-log
/// space (`-ln p`) a decaying exponential tail is a straight line, so the tail needs only a
/// single slope:
///
/// ```text
/// -ln p(t) = -ln p_edge + slope * (t - t_edge)
/// ```
///
/// equivalently `p(t) = p_edge * exp(-slope * (t - t_edge))` in plain-probability space, where
/// `p_edge = exp(-y_edge)` is recovered from the grid's own stored edge ordinate `y_edge` and
/// `t_edge` is its coordinate.
///
/// The law stores only `slope` and re-reads the current edge on evaluation. This is deliberate:
/// a soft edge is *movable* -- re-windowing and resampling shift it -- so anchoring the law to
/// the live grid edge keeps it valid across regridding, whereas a stored absolute anchor would
/// go stale. This is the essential difference from
/// [`HardApproachLaw`](crate::HardApproachLaw), whose boundary is an immovable physical fact and
/// carries an absolute anchor.
///
/// `slope` is signed so the tail decays away from support: negative on a left tail, positive on
/// a right tail. A flat tail is `slope == 0`, reached only by clamping a degenerate fit; it is
/// non-integrable, while a genuine decaying tail has finite mass.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SoftTailLaw {
  /// Signed neg-log slope `d(-ln p)/dt`: negative on a left tail, positive on a right tail.
  pub slope: f64,
}

impl SoftTailLaw {
  /// Fit a soft-tail slope from the outermost grid points adjacent to an edge.
  ///
  /// Least-squares regression of `(t, y)` over the `n_fit` points nearest the edge, where `y` is
  /// the stored neg-log ordinate, gives the neg-log slope directly. The slope is clamped so the
  /// tail decays (never grows) away from support: `slope <= 0` on the left, `slope >= 0` on the
  /// right. A wrong-sign fit means the outermost points trend back toward the peak; clamping to
  /// `0` yields a flat tail rather than one that manufactures probability outward. This reproduces
  /// v0's tail guard (`node_interpolator.py`).
  ///
  /// The fit reads stored ordinates directly, so it is valid under `NegLog` storage with no
  /// conversion, exactly like [`HardApproachLaw::fit_log_power_law`](crate::HardApproachLaw::fit_log_power_law).
  ///
  /// Returns an error when fewer than two finite points are available near the edge.
  pub fn fit(grid_fn: &GridFn<f64>, side: Side, n_fit: usize) -> Result<Self, Report> {
    let n = grid_fn.n_points();
    let n_fit = n_fit.min(n);

    let t_edge = match side {
      Side::Left => grid_fn.grid().x_at(0),
      Side::Right => grid_fn.grid().x_at(n - 1),
    };

    // Regress against `t - t_edge` rather than `t`: the slope is shift-invariant, and centering on
    // the edge keeps the normal-equation terms well conditioned when `t_edge` is large.
    let (ts, ys): (Vec<f64>, Vec<f64>) = (0..n_fit)
      .map(|i| match side {
        Side::Left => i,
        Side::Right => n - 1 - i,
      })
      .filter_map(|idx| {
        let y = grid_fn.y()[idx];
        y.is_finite().then(|| (grid_fn.grid().x_at(idx) - t_edge, y))
      })
      .collect();

    if ts.len() < 2 {
      return make_error!(
        "Soft-tail fit on the {side:?} side needs at least two finite grid points near the edge, found {}",
        ts.len()
      );
    }

    let slope_raw = LineFit::least_squares(&ts, &ys).slope;

    // Decay guard: keep only a slope that decays away from support.
    let slope = match side {
      Side::Left => slope_raw.min(0.0),
      Side::Right => slope_raw.max(0.0),
    };

    Ok(SoftTailLaw { slope })
  }

  /// Evaluate the tail in neg-log at `t`, anchored on the live grid `edge`.
  ///
  /// `edge` is the grid's current edge (its coordinate `edge.t` and stored neg-log ordinate
  /// `edge.y`). Returns `edge.y + slope * (t - edge.t)`, the straight line that a decaying
  /// exponential tail becomes in neg-log space, so the tail meets the grid continuously at the edge.
  pub fn eval(&self, edge: GridEdge, t: f64) -> f64 {
    edge.y + self.slope * (t - edge.t)
  }

  /// Tail mass beyond the edge, `exp(-y_edge) / |slope|`, in plain probability.
  ///
  /// Closed form of the half-line integral `∫ p_edge exp(-slope (t - t_edge)) dt`, with the edge
  /// probability recovered from its stored neg-log ordinate as `p_edge = exp(-y_edge)`. A flat
  /// tail (`slope == 0`) is non-integrable and returns `+inf`; a genuine decaying tail is finite.
  pub fn mass(&self, y_edge: f64) -> f64 {
    (-y_edge).exp() / self.slope.abs()
  }

  /// Compose two soft tails under multiplication: the neg-log slopes add.
  ///
  /// Multiplication is addition in neg-log space, so the product tail's slope is the sum of the
  /// operand slopes. Because both operands and the product are evaluated against the same live
  /// grid edge, the result needs only the summed slope -- there is no anchor to reconcile.
  #[must_use]
  pub fn compose_multiply(&self, other: &SoftTailLaw) -> SoftTailLaw {
    SoftTailLaw {
      slope: self.slope + other.slope,
    }
  }

  /// Reflect the tail for `f(t) -> f(-t)`: the slope flips sign.
  #[must_use]
  pub fn negate_arg(&self) -> Self {
    SoftTailLaw { slope: -self.slope }
  }
}
