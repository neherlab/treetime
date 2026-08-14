use crate::GridFn;
use serde::{Deserialize, Serialize};

/// Hard-boundary approach law in neg-log space.
///
/// Evaluates the neg-log density in the sub-grid gap `[t_hard, t_first)` between a hard boundary
/// and the nearest grid point:
///
/// ```text
/// y(t) = a - b * ln|t - t_hard| + slope * (t - t_hard)
/// ```
///
/// Two regimes arise from the branch-length Gamma likelihood `p(t) ~ t^n * exp(-mu*t)`:
///
/// - `n >= 1` mutations: `y(t) = -n * ln(t) + mu*t + C`, so `b = n > 0` and `slope = mu`.
///   The density vanishes at the boundary (`y -> +inf`). The log term dominates, and the linear
///   term contributes a small correction over the sub-mutation gap. The fit sets `slope = 0`
///   because the gap is too narrow for the linear term to be resolved.
///
/// - `n = 0` mutations: `y(t) = mu*t + C`, so `b = 0` and `slope = mu`.
///   The density is maximal at the boundary (`y` is finite). The fit recovers the exact linear
///   form, preserving the mode on the boundary.
///
/// Under multiplication (addition in neg-log) all three parameters add:
///
/// ```text
/// (a1 - b1*ln|dt| + s1*dt) + (a2 - b2*ln|dt| + s2*dt) = (a1+a2) - (b1+b2)*ln|dt| + (s1+s2)*dt
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HardApproachLaw {
  /// Location of the hard boundary (e.g. t=0 for branch lengths).
  pub t_hard: f64,
  /// Neg-log anchor: the `a` in `y(t) = a - b*ln|dt| + slope*dt`.
  /// For `b = 0` this is the neg-log ordinate at the boundary itself.
  /// For `b > 0` it is the intercept of the log-distance regression projected to `ln|dt| = 0`.
  pub a: f64,
  /// Power-law exponent `b >= 0`. Equals the mutation count `n` for a Gamma branch-length
  /// density. `b = 0` means the density is finite at the boundary; `b > 0` means it vanishes.
  pub b: f64,
  /// Linear slope `d(-ln p)/dt` near the boundary. For the initial fit this is set from a
  /// linear regression when `b = 0`, and set to `0` when `b > 0` (the log term dominates the
  /// narrow sub-grid gap). Composition adds slopes, so a product of zero- and nonzero-mutation
  /// messages carries both `b > 0` and `slope > 0`.
  pub slope: f64,
}

impl HardApproachLaw {
  /// Fit a hard-boundary approach law from the innermost grid points.
  ///
  /// The fit reads stored ordinates directly (neg-log values under `NegLog` storage). Two-stage:
  ///
  /// 1. Log-distance regression on `(ln|t - t_hard|, y_stored)` gives `-b` as the slope. If
  ///    `b > threshold`, the power-law form applies: `a` is the intercept projected to
  ///    `ln|dt| = 0`, `slope` is set to `0`.
  ///
  /// 2. If `b <= threshold` (clamped to 0), the density is finite at the boundary: refit
  ///    linearly on `(t, y_stored)` to recover `(slope, a)` exactly.
  ///
  /// `side` selects which end of the grid is near the hard boundary:
  /// - `Side::Left`: the hard boundary is to the left of `x_min`, fit from the leftmost points
  /// - `Side::Right`: the hard boundary is to the right of `x_max`, fit from the rightmost points
  ///
  /// Returns `None` when fewer than two finite points are available near the boundary, or when
  /// the fit is not finite.
  pub fn fit(grid_fn: &GridFn<f64>, t_hard: f64, side: Side, n_fit: usize) -> Option<Self> {
    let n = grid_fn.n_points();
    let n_fit = n_fit.min(n);

    let points: Vec<(f64, f64)> = (0..n_fit)
      .map(|i| match side {
        Side::Left => i,
        Side::Right => n - 1 - i,
      })
      .filter_map(|idx| {
        let y = grid_fn.y()[idx];
        y.is_finite().then(|| (grid_fn.grid().x_at(idx), y))
      })
      .collect();

    if points.len() < 2 {
      return None;
    }

    // Stage 1: log-distance regression to detect power-law exponent.
    // Transform: (ln|t - t_hard|, y_stored). The neg-log power law y = a - b*ln|dt| is linear
    // in ln|dt| with slope = -b.
    let log_dist_points: Vec<(f64, f64)> = points
      .iter()
      .filter_map(|&(t, y)| {
        let dt = (t - t_hard).abs();
        (dt > 0.0).then(|| (dt.ln(), y))
      })
      .collect();

    if log_dist_points.len() < 2 {
      return None;
    }

    let (xs, ys): (Vec<f64>, Vec<f64>) = log_dist_points.into_iter().unzip();
    let (neg_b_raw, _) = least_squares_fit(&xs, &ys);
    let b = (-neg_b_raw).max(0.0);

    const B_THRESHOLD: f64 = 0.01;

    if b > B_THRESHOLD {
      // Power-law case (n >= 1 mutations). Re-derive the intercept conditioned on the clamped
      // exponent: `a = mean(y) + b * mean(ln|dt|)`.
      let x_mean = xs.iter().sum::<f64>() / xs.len() as f64;
      let y_mean = ys.iter().sum::<f64>() / ys.len() as f64;
      let a = y_mean + b * x_mean;

      if !a.is_finite() {
        return None;
      }

      Some(HardApproachLaw {
        t_hard,
        a,
        b,
        slope: 0.0,
      })
    } else {
      // Linear case (n = 0 mutations). Refit on (t, y_stored) to recover the exact linear form.
      let (ts, ys): (Vec<f64>, Vec<f64>) = points.into_iter().unzip();
      let (slope, intercept) = least_squares_fit(&ts, &ys);
      let a = slope * t_hard + intercept;

      if !slope.is_finite() || !a.is_finite() {
        return None;
      }

      Some(HardApproachLaw {
        t_hard,
        a,
        b: 0.0,
        slope,
      })
    }
  }

  /// Evaluate the approach law at a point between `t_hard` and the nearest grid point.
  ///
  /// Returns the neg-log ordinate `a - b*ln|t - t_hard| + slope*(t - t_hard)`.
  /// When `b > 0` and `t == t_hard`, the density is zero (neg-log is `+inf`).
  /// When `b == 0`, the value at the boundary is `a` (finite, the mode).
  pub fn eval(&self, t: f64) -> f64 {
    let dt = (t - self.t_hard).abs();
    if dt == 0.0 {
      if self.b == 0.0 {
        return self.a;
      }
      return f64::INFINITY;
    }
    self.a - self.b * dt.ln() + self.slope * dt
  }

  /// Compose two approach laws under multiplication.
  ///
  /// Multiplication is addition in neg-log, so all three parameters add.
  /// Both laws must share the same `t_hard`.
  #[must_use]
  pub fn compose_multiply(&self, other: &HardApproachLaw) -> HardApproachLaw {
    HardApproachLaw {
      t_hard: self.t_hard,
      a: self.a + other.a,
      b: self.b + other.b,
      slope: self.slope + other.slope,
    }
  }

  /// Probability mass in the approach region between `t_hard` and `t_grid`.
  ///
  /// Closed-form integral of `exp(-y(t))` over the gap `[t_hard, t_grid]`.
  ///
  /// - `b > 0, slope = 0`: `exp(-a) * dt^(b+1) / (b+1)` (power-law in plain space)
  /// - `b = 0, slope != 0`: `exp(-a) * (1 - exp(-slope*dt)) / slope`
  /// - `b = 0, slope = 0`: `exp(-a) * dt`
  /// - `b > 0, slope > 0`: lower incomplete gamma; approximated by the `slope = 0` formula
  ///   since `slope*dt` is small over the sub-grid gap.
  pub fn mass(&self, t_grid: f64) -> f64 {
    let dt = (t_grid - self.t_hard).abs();
    let scale = (-self.a).exp();

    if self.b > 0.0 {
      // Power-law case. When slope is also nonzero (from composition), the integral is an
      // incomplete gamma function. Over the sub-grid gap slope*dt is small, so the leading
      // power-law term dominates.
      scale * dt.powf(self.b + 1.0) / (self.b + 1.0)
    } else if self.slope.abs() > 0.0 {
      // Linear case: integral of exp(-a - slope*dt') dt' from 0 to dt.
      scale * (1.0 - (-self.slope * dt).exp()) / self.slope
    } else {
      // Flat case: constant density.
      scale * dt
    }
  }

  /// Transform the approach law when the argument is negated: `f(x) -> f(-x)`.
  ///
  /// The boundary moves to `-t_hard`, the anchor `a` and exponent `b` are unchanged,
  /// and the linear slope flips sign.
  #[must_use]
  pub fn negate_arg(&self) -> Self {
    HardApproachLaw {
      t_hard: -self.t_hard,
      a: self.a,
      b: self.b,
      slope: -self.slope,
    }
  }

  /// Shift the anchor when all ordinates are offset by a constant.
  ///
  /// Under NegLog normalization, subtracting a constant `c` from every stored ordinate
  /// shifts `a` to `a - c`. The exponent `b` and slope are invariant (they depend on
  /// shape, not scale).
  #[must_use]
  pub fn shift_anchor(&self, offset: f64) -> Self {
    HardApproachLaw {
      a: self.a + offset,
      ..*self
    }
  }
}

/// Which side of the grid a hard boundary is on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
  Left,
  Right,
}

/// Simple least-squares linear regression: y = slope * x + intercept.
pub(crate) fn least_squares_fit(xs: &[f64], ys: &[f64]) -> (f64, f64) {
  let n = xs.len() as f64;
  let sum_x: f64 = xs.iter().sum();
  let sum_y: f64 = ys.iter().sum();
  let sum_xx: f64 = xs.iter().map(|x| x * x).sum();
  let sum_xy: f64 = xs.iter().zip(ys).map(|(x, y)| x * y).sum();

  let sum_x_sq = sum_x * sum_x;
  let denom = n * sum_xx - sum_x_sq;
  if denom.abs() < 1e-30 {
    return (0.0, sum_y / n);
  }

  let slope = (n * sum_xy - sum_x * sum_y) / denom;
  let intercept = (sum_y - slope * sum_x) / n;
  (slope, intercept)
}
