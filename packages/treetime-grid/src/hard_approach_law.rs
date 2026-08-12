use crate::GridFn;
use serde::{Deserialize, Serialize};

/// Power-law interpolation between a hard boundary and the nearest grid point.
///
/// Near a hard boundary the density follows `p(t) = C * |t - t_hard|^b` in
/// plain-probability space, equivalently `y(t) = a - b * log|t - t_hard|` in
/// neg-log space where `C = exp(-a)`.
///
/// - `b = 0`: constant at the boundary (zero-mutation branch, density is maximal
///   at t=0). The grid endpoint stores the finite boundary value directly.
/// - `b > 0`: power-law decay toward zero (n-mutation branch, `p(t) ~ t^n`). The
///   grid stays strictly inside the hard boundary and the approach law covers the
///   gap `[t_hard, t_first)`.
///
/// Under multiplication the exponents add (`b_result = b_a + b_b`) because the
/// product of two power laws is a power law with summed exponents.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HardApproachLaw {
  /// Location of the hard boundary (e.g. t=0 for branch lengths).
  pub t_hard: f64,
  /// Amplitude coefficient C in `p(t) = C * |t - t_hard|^b`. Always positive.
  pub coeff: f64,
  /// Power-law exponent b >= 0.
  pub exponent: f64,
}

impl HardApproachLaw {
  /// Fit an approach law from the innermost grid points adjacent to a hard boundary.
  ///
  /// Uses least-squares regression on `(log|t - t_hard|, log(p(t)))` from the `n_fit`
  /// innermost points. The slope gives `b`; the coefficient `C` is `exp(a)`, where the
  /// intercept `a` is re-derived from the (possibly clamped) exponent so that `C` stays
  /// consistent with `b`.
  ///
  /// Precondition: `grid_fn.y()` holds *plain-probability* values. The fit takes their
  /// natural log, so it is valid only under the `Plain` y-axis policy. Under `NegLog`
  /// storage the caller must convert to plain probability first (`-Y::to_neg_log(y)` yields
  /// the log-probability directly).
  ///
  /// `side` selects which end of the grid is near the hard boundary:
  /// - `Side::Left`: the hard boundary is to the left of `x_min`, fit from the leftmost points
  /// - `Side::Right`: the hard boundary is to the right of `x_max`, fit from the rightmost points
  ///
  /// Returns `None` when the fit is not meaningful (all y values zero, or fewer than 2
  /// positive points near the boundary).
  pub fn fit(grid_fn: &GridFn<f64>, t_hard: f64, side: Side, n_fit: usize) -> Option<Self> {
    let n = grid_fn.n_points();
    let n_fit = n_fit.min(n);

    // NOTE: `grid_fn.y()[i]` is read as a plain-probability value and fed to `.ln()`.
    // This is correct only while the source distribution uses the `Plain` y-axis policy. When
    // timetree distributions switch to `NegLog`, the stored values become neg-log and this
    // `.ln()` silently produces a wrong law; the caller must pass plain-probability values then.
    let (xs, ys): (Vec<f64>, Vec<f64>) = match side {
      Side::Left => (0..n_fit)
        .filter_map(|i| {
          let x = grid_fn.grid().x_at(i);
          let y = grid_fn.y()[i];
          (y > 0.0).then(|| ((x - t_hard).abs().ln(), y.ln()))
        })
        .collect(),
      Side::Right => (0..n_fit)
        .filter_map(|i| {
          let idx = n - 1 - i;
          let x = grid_fn.grid().x_at(idx);
          let y = grid_fn.y()[idx];
          (y > 0.0).then(|| ((t_hard - x).abs().ln(), y.ln()))
        })
        .collect(),
    };

    if xs.len() < 2 {
      return None;
    }

    let (slope, _joint_intercept) = least_squares_fit(&xs, &ys);

    // b >= 0: a negative slope means the density increases into the boundary faster
    // than any power law, which is unphysical; clamp to b = 0.
    let exponent = slope.max(0.0);

    // NOTE: re-derive the intercept conditioned on the clamped exponent instead of
    // reusing the joint OLS intercept. The conditional-least-squares intercept for a fixed
    // slope `b` is `mean(y) - b * mean(x)`. Unclamped (`b == slope`) it equals the joint
    // intercept (no change); clamped (`b == 0`) it collapses to `mean(log p)`, the correct
    // horizontal fit. Reusing the joint intercept after clamping biased `coeff` low by
    // `exp(-slope * mean(x))`.
    let x_mean = xs.iter().sum::<f64>() / xs.len() as f64;
    let y_mean = ys.iter().sum::<f64>() / ys.len() as f64;
    let coeff = (y_mean - exponent * x_mean).exp();

    if !coeff.is_finite() || coeff <= 0.0 {
      return None;
    }

    Some(HardApproachLaw {
      t_hard,
      coeff,
      exponent,
    })
  }

  /// Evaluate the approach law at a point between `t_hard` and the nearest grid point.
  ///
  /// Returns `C * |t - t_hard|^b` in plain-probability space.
  pub fn eval(&self, t: f64) -> f64 {
    let dt = (t - self.t_hard).abs();
    if dt == 0.0 {
      if self.exponent == 0.0 {
        return self.coeff;
      }
      return 0.0;
    }
    self.coeff * dt.powf(self.exponent)
  }

  /// Compose two approach laws under multiplication (exponents add).
  ///
  /// The product `p_a(t) * p_b(t) = C_a * dt^b_a * C_b * dt^b_b = (C_a * C_b) * dt^(b_a + b_b)`.
  /// Both laws must share the same `t_hard`.
  #[must_use]
  pub fn compose_multiply(&self, other: &HardApproachLaw) -> HardApproachLaw {
    HardApproachLaw {
      t_hard: self.t_hard,
      coeff: self.coeff * other.coeff,
      exponent: self.exponent + other.exponent,
    }
  }

  /// Mass in the approach region between `t_hard` and `t_grid` (the nearest grid point).
  ///
  /// Closed form: `C * |t_grid - t_hard|^(b+1) / (b + 1)`.
  pub fn mass(&self, t_grid: f64) -> f64 {
    let dt = (t_grid - self.t_hard).abs();
    self.coeff * dt.powf(self.exponent + 1.0) / (self.exponent + 1.0)
  }

  /// Transform the approach law when the argument is negated: `f(x) -> f(-x)`.
  /// The hard boundary moves to `-t_hard`.
  #[must_use]
  pub fn negate_arg(&self) -> Self {
    HardApproachLaw {
      t_hard: -self.t_hard,
      coeff: self.coeff,
      exponent: self.exponent,
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
