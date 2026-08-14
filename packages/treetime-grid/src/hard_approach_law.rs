use crate::GridFn;
use serde::{Deserialize, Serialize};

/// Log-distance regression slope above which the boundary is treated as divergent (power-law,
/// `b > 0`, `slope = 0`); at or below it the branch is finite (`b = 0`) and the slope is refit.
const B_THRESHOLD: f64 = 0.01;

/// Approach law for a *hard* grid boundary, in negative-log space.
///
/// Between a hard boundary at `t_hard` and the nearest grid point, the branch-length density
/// follows the Gamma likelihood `p(t) ~ |t - t_hard|^b * exp(-slope*|t - t_hard|)`, so in
/// negative-log space it is the two-term law
///
/// ```text
/// y(t) = y_edge - b*ln( |t - t_hard| / |t_edge - t_hard| ) + slope*(t - t_edge)
/// ```
///
/// where `y_edge` is the grid's stored neg-log edge ordinate (`y[0]` on the left, `y[n-1]` on the
/// right) and `t_edge` its coordinate. The law is *edge-relative*: it stores only the shape
/// parameters `b` and `slope` and reads the live grid edge on evaluation, exactly like
/// [`SoftTailLaw`](crate::SoftTailLaw). Reading the live edge keeps the law valid across
/// re-windowing and resampling with no refit and no absolute anchor to keep in sync: a vertical
/// shift of every ordinate (peak-normalization) is absorbed through `y_edge`, and a scale
/// (`p -> p^factor`) scales both shape parameters.
///
/// Two regimes arise from the branch-length likelihood `p(t) ~ t^n * exp(-mu*t)` (`n` mutations):
///
/// - `n >= 1` (divergent): `y(t) = -n*ln(t) + mu*t + C`, so `b = n > 0`. The density vanishes at the
///   boundary (`y -> +inf`). The log term dominates; the linear term is a small correction over the
///   sub-mutation gap, so the fit sets `slope = 0`.
/// - `n = 0` (finite): `y(t) = mu*t + C`, so `b = 0` and `slope = mu`. The density is finite and
///   maximal at the boundary; the linear term carries the mode, which sits exactly on the edge.
///
/// A single law covers both: `n = 0` is not a special case, it is `b = 0`, where the log term
/// vanishes and the law is the line `y_edge + slope*(t - t_edge)`. This matches the derived Part C
/// truth `y = -n*ln(t) + mu*t + C`, of which the flat single-term proposed form (`a - b*ln|dt|`) is
/// the `slope = 0` restriction that cannot represent the `n = 0` mode.
///
/// The boundary location `t_hard` is an immovable physical fact (e.g. `t = 0` for branch lengths)
/// and is stored absolutely; `b` and `slope` are the shape.
///
/// Under multiplication (addition in neg-log) both shape parameters add:
/// `b_result = b_a + b_b`, `slope_result = slope_a + slope_b`. A product of a divergent and a
/// finite message therefore carries both a `b > 0` and a `slope > 0` term, which is why the two
/// terms must live in one law rather than in separate variants that would not close under
/// multiplication.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HardApproachLaw {
  /// Location of the hard boundary (e.g. `t = 0` for branch lengths).
  pub t_hard: f64,
  /// Power-law exponent `b >= 0`. Equals the mutation count `n` for a Gamma branch-length density.
  /// `b > 0` means the density vanishes at the boundary (`y -> +inf`); `b = 0` means it is finite.
  pub b: f64,
  /// Linear slope `d(-ln p)/dt` near the boundary. Carries the `n = 0` mode (`slope = mu`). The
  /// initial fit sets it from a linear regression when `b = 0`, and to `0` when `b > 0` (the log
  /// term dominates the narrow sub-grid gap). Composition adds slopes, so a product of zero- and
  /// nonzero-mutation messages carries both `b > 0` and `slope > 0`.
  pub slope: f64,
}

impl HardApproachLaw {
  /// Fit a hard-boundary approach law from the innermost grid points.
  ///
  /// The fit reads stored ordinates directly (neg-log values under `NegLog` storage). Two-stage:
  ///
  /// 1. Log-distance regression on `(ln|t - t_hard|, y)` gives `-b` as the slope. If `b >
  ///    threshold`, the power-law form applies and `slope` is set to `0` (the log term dominates the
  ///    narrow sub-grid gap).
  /// 2. If `b <= threshold` (clamped to `0`), the density is finite at the boundary: refit linearly
  ///    on `(t, y)` to recover the exact `slope`, which carries the mode on the boundary.
  ///
  /// The intercept is discarded in both stages because the law is edge-relative; only the shape
  /// parameters `b` and `slope` are retained. `b` is clamped to `b >= 0` (a negative fit means the
  /// density grows into the boundary faster than any power law, which is unphysical here).
  ///
  /// `side` selects which end of the grid is near the hard boundary:
  /// - `Side::Left`: the hard boundary is to the left of `x_min`, fit from the leftmost points
  /// - `Side::Right`: the hard boundary is to the right of `x_max`, fit from the rightmost points
  ///
  /// Returns `None` when fewer than two finite points are available near the boundary, or when the
  /// fit is not finite.
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

    // Stage 1: log-distance regression to detect the power-law exponent. Transform to
    // (ln|t - t_hard|, y); the neg-log power law y = a - b*ln|dt| is linear in ln|dt| with slope
    // = -b. The intercept is discarded (edge-relative).
    let (xs, ys): (Vec<f64>, Vec<f64>) = points
      .iter()
      .filter_map(|&(t, y)| {
        let dt = (t - t_hard).abs();
        (dt > 0.0).then(|| (dt.ln(), y))
      })
      .unzip();

    if xs.len() < 2 {
      return None;
    }

    let (neg_b_raw, _) = least_squares_fit(&xs, &ys);
    let b = (-neg_b_raw).max(0.0);

    if !b.is_finite() {
      return None;
    }

    if b > B_THRESHOLD {
      // Power-law case (n >= 1 mutations): the log term dominates the narrow gap, so drop the
      // linear correction.
      Some(HardApproachLaw { t_hard, b, slope: 0.0 })
    } else {
      // Linear case (n = 0 mutations): refit on (t, y) to recover the exact slope, which carries
      // the mode on the boundary.
      let (ts, ys): (Vec<f64>, Vec<f64>) = points.into_iter().unzip();
      let (slope, _intercept) = least_squares_fit(&ts, &ys);

      if !slope.is_finite() {
        return None;
      }

      Some(HardApproachLaw { t_hard, b: 0.0, slope })
    }
  }

  /// Evaluate the approach law in neg-log at `t`, anchored on the live grid edge.
  ///
  /// `y_edge` is the grid's stored neg-log edge ordinate and `t_edge` its coordinate. Returns
  /// `y_edge - b*ln(|t - t_hard| / |t_edge - t_hard|) + slope*(t - t_edge)`, so the law meets the
  /// grid continuously at the edge. At `t == t_hard` a divergent law (`b > 0`) returns `+inf`
  /// (density zero at the boundary); a finite law (`b = 0`) returns the line's value there,
  /// `y_edge + slope*(t_hard - t_edge)`.
  pub fn eval(&self, y_edge: f64, t_edge: f64, t: f64) -> f64 {
    let dt = (t - self.t_hard).abs();
    let dt_edge = (t_edge - self.t_hard).abs();
    let log_term = if self.b == 0.0 {
      0.0
    } else if dt == 0.0 {
      return f64::INFINITY;
    } else {
      self.b * (dt / dt_edge).ln()
    };
    y_edge - log_term + self.slope * (t - t_edge)
  }

  /// Probability mass in the approach region between `t_hard` and the grid edge, in plain
  /// probability.
  ///
  /// Closed form of `integral_{t_hard}^{t_edge} exp(-y(t)) dt` with the edge probability recovered
  /// from its stored neg-log ordinate as `p_edge = exp(-y_edge)`:
  ///
  /// - `b > 0`: `p_edge * |t_edge - t_hard| / (b + 1)` (power law; the `slope` correction is
  ///   negligible over the sub-grid gap and is dropped, mirroring the fit).
  /// - `b = 0, slope != 0`: `p_edge * (exp(slope*dt_edge) - 1) / slope`.
  /// - `b = 0, slope = 0`: `p_edge * dt_edge`.
  pub fn mass(&self, y_edge: f64, t_edge: f64) -> f64 {
    let dt_edge = (t_edge - self.t_hard).abs();
    let p_edge = (-y_edge).exp();

    if self.b > 0.0 {
      p_edge * dt_edge / (self.b + 1.0)
    } else if self.slope != 0.0 {
      p_edge * (self.slope * dt_edge).exp_m1() / self.slope
    } else {
      p_edge * dt_edge
    }
  }

  /// Compose two approach laws under multiplication: the shape parameters add.
  ///
  /// Multiplication is addition in neg-log space, so the product law's exponent and slope are the
  /// sums of the operand values. Both laws must share the same `t_hard`, and both are evaluated
  /// against the same live grid edge, so there is no anchor to reconcile.
  #[must_use]
  pub fn compose_multiply(&self, other: &HardApproachLaw) -> HardApproachLaw {
    HardApproachLaw {
      t_hard: self.t_hard,
      b: self.b + other.b,
      slope: self.slope + other.slope,
    }
  }

  /// Transform the approach law when the argument is negated: `f(x) -> f(-x)`.
  ///
  /// The boundary moves to `-t_hard`; the exponent `b` is unchanged and the linear slope flips sign.
  #[must_use]
  pub fn negate_arg(&self) -> Self {
    HardApproachLaw {
      t_hard: -self.t_hard,
      b: self.b,
      slope: -self.slope,
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
