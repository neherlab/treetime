use crate::{GridEdge, GridFn};
use eyre::Report;
use serde::{Deserialize, Serialize};
use treetime_utils::least_squares::LineFit;
use treetime_utils::make_error;

/// Approach law for a *hard* grid boundary, in negative-log space.
///
/// A hard boundary at `t_hard` terminates the domain: probability is zero beyond it. Between the
/// boundary and the nearest grid point the density follows a one-sided power law
/// `p(t) ~ (t - t_hard)^b` with a single exponent `b >= 0`. For a branch-length factor
/// `L(t) ~ (mu*t)^n exp(-mu*t)` the leading behaviour as `t -> t_hard` is exactly this power law with
/// `b = n` (mutation count), the `exp(-mu*t)` factor contributing only a sub-leading term that
/// vanishes at the boundary. In negative-log space the law is
///
/// ```text
/// y(t) = a - b * ln|t - t_hard|
/// ```
///
/// a straight line in `ln|t - t_hard|`. The exponent `b` is the only shape parameter; the intercept
/// `a` is not stored because the law is *edge-relative* (see below).
///
/// `b > 0` is the divergent case (a branch with mutations, or an indel): the density vanishes at the
/// boundary and the neg-log diverges to `+inf`. `b = 0` is the flat case (`y = a`), where the density
/// is finite at the boundary. A finite mode sitting exactly on a hard boundary -- the common
/// zero-mutation branch, whose density is maximal at `t = 0` -- is not represented by this law at all:
/// the producer places the boundary as an exact grid endpoint ([`BoundaryBehavior::Hard`]), so the
/// grid carries the finite ordinate and the mode directly, with no sub-grid gap to approximate.
/// `b = 0` here is only the continuum value the fit returns when the innermost points are flat.
///
/// The law is *edge-relative*: it stores only `t_hard` and `b`, and reads the live grid edge
/// ([`GridEdge`](crate::GridEdge)) on evaluation, exactly like [`SoftTailLaw`](crate::SoftTailLaw).
/// The boundary location `t_hard` is an immovable physical fact (e.g. `t = 0` for branch lengths); the
/// ordinate is read live, so a peak-normalization shift, a resample, or a re-window needs no refit.
///
/// [`BoundaryBehavior::Hard`]: crate::BoundaryBehavior::Hard
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HardApproachLaw {
  /// Location of the hard boundary (e.g. `t = 0` for branch lengths).
  pub t_hard: f64,
  /// Power-law exponent `b >= 0` of the boundary approach `p(t) ~ (t - t_hard)^b`. `b > 0` diverges in
  /// neg-log at the boundary; `b = 0` is flat (finite boundary density).
  pub b: f64,
}

impl HardApproachLaw {
  /// Fit the power-law exponent `b >= 0` from the innermost `n_fit` grid points nearest the boundary.
  ///
  /// In neg-log the approach is `y = a - b*ln|t - t_hard|`, linear in `ln|t - t_hard|` with slope
  /// `-b`. One least-squares regression on `(ln|t - t_hard|, y)` over the innermost finite points off
  /// the boundary recovers the slope; the intercept `a` is discarded (edge-relative). A wrong-sign fit
  /// (`b < 0`: density increasing into the boundary faster than any power law, unphysical here) clamps
  /// to `b = 0`, a flat shape, rather than fabricating a boundary singularity the data does not
  /// support.
  ///
  /// `Side` selects the grid end nearest the boundary. Returns an error on fewer than two finite
  /// innermost points off the boundary, or a non-finite result.
  pub fn fit(grid_fn: &GridFn<f64>, t_hard: f64, side: Side, n_fit: usize) -> Result<Self, Report> {
    let n = grid_fn.n_points();
    let n_fit = n_fit.min(n);
    let (xs, ys): (Vec<f64>, Vec<f64>) = (0..n_fit)
      .map(|i| match side {
        Side::Left => i,
        Side::Right => n - 1 - i,
      })
      .filter_map(|idx| {
        let y = grid_fn.y()[idx];
        let dt = (grid_fn.grid().x_at(idx) - t_hard).abs();
        (y.is_finite() && dt > 0.0).then(|| (dt.ln(), y))
      })
      .collect();

    if xs.len() < 2 {
      return make_error!(
        "Hard-boundary power-law fit on the {side:?} side needs at least two finite grid points off \
         the boundary t_hard={t_hard}, found {}",
        xs.len()
      );
    }

    let neg_b_raw = LineFit::least_squares(&xs, &ys).slope;
    let b = (-neg_b_raw).max(0.0);
    if !b.is_finite() {
      return make_error!("Hard-boundary power-law fit on the {side:?} side produced a non-finite exponent");
    }
    Ok(HardApproachLaw { t_hard, b })
  }

  /// Evaluate the approach law in neg-log at `t`, anchored on the live grid `edge`.
  ///
  /// `edge` is the grid's current edge (its coordinate `edge.t` and stored neg-log ordinate
  /// `edge.y`). The law is `edge.y - b*ln(|t - t_hard| / |edge.t - t_hard|)`, so it meets the grid
  /// continuously at the edge (`t == edge.t` gives `edge.y`). At the boundary `t == t_hard` it is
  /// `+inf` for `b > 0` (the density is zero there) and the finite `edge.y` for `b == 0` (flat).
  pub fn eval(&self, edge: GridEdge, t: f64) -> f64 {
    let dt = (t - self.t_hard).abs();
    if dt == 0.0 {
      return if self.b == 0.0 { edge.y } else { f64::INFINITY };
    }
    let dt_edge = (edge.t - self.t_hard).abs();
    edge.y - self.b * (dt / dt_edge).ln()
  }

  /// Probability mass in the approach region between `t_hard` and the grid edge, in plain probability.
  ///
  /// Closed form of `integral_{t_hard}^{edge.t} exp(-y(t)) dt` with the edge probability recovered from
  /// the anchor's stored neg-log ordinate as `p_edge = exp(-edge.y)`:
  ///
  /// ```text
  /// p_edge * |edge.t - t_hard| / (b + 1)
  /// ```
  ///
  /// which reduces to `p_edge * dt_edge` in the flat case `b = 0`.
  pub fn mass(&self, edge: GridEdge) -> f64 {
    let dt_edge = (edge.t - self.t_hard).abs();
    let p_edge = (-edge.y).exp();
    p_edge * dt_edge / (self.b + 1.0)
  }

  /// Scale the approach law when every neg-log ordinate is multiplied by `factor` (`p -> p^factor`).
  ///
  /// The exponent scales by `factor`; the boundary location is unchanged.
  #[must_use]
  pub fn scale(&self, factor: f64) -> HardApproachLaw {
    HardApproachLaw {
      t_hard: self.t_hard,
      b: self.b * factor,
    }
  }

  /// Transform the approach law when the argument is negated: `f(x) -> f(-x)`.
  ///
  /// The boundary moves to `-t_hard`; the exponent `b` is unchanged.
  #[must_use]
  pub fn negate_arg(&self) -> Self {
    HardApproachLaw {
      t_hard: -self.t_hard,
      b: self.b,
    }
  }
}

/// Which side of the grid a hard boundary is on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
  Left,
  Right,
}
