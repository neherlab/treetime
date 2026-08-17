use crate::GridFn;
use serde::{Deserialize, Serialize};

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
  /// Power-law exponent `b >= 0`. `b > 0`: density vanishes at the boundary (`y -> +inf`), a power
  /// law fitted from the grid; a fitted continuous approximation of the mutation count, not an
  /// integer. `b = 0`: finite boundary, the law is a straight line.
  pub b: f64,
  /// Linear slope of the neg-log law near the boundary. Finite boundary: the slope of the line from
  /// the boundary to the grid edge. Divergent boundary: `0`. Composition adds slopes.
  pub slope: f64,
}

impl HardApproachLaw {
  /// Fit the straight-line approach for a finite boundary density (`n = 0`) over `[t_hard, t_edge)`.
  ///
  /// A zero-mutation branch has a finite, maximal density at the boundary, so the neg-log is a
  /// straight line. The law is that exact line through the boundary point `(t_hard, y_hard)` and the
  /// grid edge -- two points, no regression: `b = 0`, `slope = (y_edge - y_hard) / (t_edge -
  /// t_hard)`. `y_hard = -ln p(t_hard)` is the boundary neg-log ordinate, supplied by the producer
  /// (the grid samples strictly inside `min_bl > 0`, so no grid point sits at `t_hard`).
  ///
  /// `Side` selects the grid end nearest the boundary. `None` on an empty grid, a non-finite edge
  /// ordinate, or a non-finite slope.
  pub fn fit_linear(grid_fn: &GridFn<f64>, t_hard: f64, side: Side, y_hard: f64) -> Option<Self> {
    let (t_edge, y_edge) = edge_ordinate(grid_fn, side)?;
    let slope = (y_edge - y_hard) / (t_edge - t_hard);
    slope.is_finite().then_some(HardApproachLaw { t_hard, b: 0.0, slope })
  }

  /// Fit the power-law approach for a divergent boundary density (`n >= 1`, or an indel) over
  /// `[t_hard, t_edge)`.
  ///
  /// The density vanishes at the boundary, so `-ln p` diverges. Fit `b >= 0` from the innermost
  /// `n_fit` points by one log-distance regression on `(ln|t - t_hard|, y)`: `y = a - b*ln|dt|` is
  /// linear in `ln|dt|` with slope `-b`; `slope = 0`. The intercept is discarded (edge-relative).
  ///
  /// `Side` selects the grid end nearest the boundary. `None` on an empty grid, a non-finite edge
  /// ordinate, fewer than two finite innermost points, or a non-finite result.
  pub fn fit_log_power_law(grid_fn: &GridFn<f64>, t_hard: f64, side: Side, n_fit: usize) -> Option<Self> {
    let n = grid_fn.n_points();
    edge_ordinate(grid_fn, side)?;

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
      return None;
    }

    let (neg_b_raw, _) = least_squares_fit(&xs, &ys);
    let b = (-neg_b_raw).max(0.0);
    b.is_finite().then_some(HardApproachLaw { t_hard, b, slope: 0.0 })
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

/// Grid edge nearest the boundary as `(t_edge, y_edge)`.
///
/// `None` on an empty grid or a non-finite edge ordinate, which cannot anchor an edge-relative law.
fn edge_ordinate(grid_fn: &GridFn<f64>, side: Side) -> Option<(f64, f64)> {
  let n = grid_fn.n_points();
  if n == 0 {
    return None;
  }
  let edge_idx = match side {
    Side::Left => 0,
    Side::Right => n - 1,
  };
  let t_edge = grid_fn.grid().x_at(edge_idx);
  let y_edge = grid_fn.y()[edge_idx];
  y_edge.is_finite().then_some((t_edge, y_edge))
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
