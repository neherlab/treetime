use crate::{GridEdge, GridFn};
use eyre::Report;
use serde::{Deserialize, Serialize};
use treetime_utils::least_squares::LineFit;
use treetime_utils::make_error;

/// Approach law for a *hard* grid boundary, in negative-log space.
///
/// Between a hard boundary at `t_hard` and the nearest grid point the neg-log density follows one of
/// two [regimes](Approach): a finite line (a zero-mutation branch, whose mode sits on the boundary)
/// or a power-law divergence (a branch with mutations, whose density vanishes at the boundary). The
/// derivation and the surrounding tail framework are recorded in
/// `kb/decisions/distribution-tails-and-arithmetic.md`.
///
/// The law is *edge-relative*: it stores only shape (`t_hard` and the regime) and reads the live
/// grid edge ([`GridEdge`](crate::GridEdge)) on evaluation, exactly like
/// [`SoftTailLaw`](crate::SoftTailLaw). The
/// boundary location `t_hard` is an immovable physical fact (e.g. `t = 0` for branch lengths); the
/// ordinate is read live, so a peak-normalization shift, a resample, or a re-window needs no refit.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HardApproachLaw {
  /// Location of the hard boundary (e.g. `t = 0` for branch lengths).
  pub t_hard: f64,
  /// The neg-log approach shape near the boundary. See [`Approach`].
  pub shape: Approach,
}

impl HardApproachLaw {
  /// Fit the straight-line approach for a finite boundary density (`n = 0`) over `[t_hard, t_edge)`.
  ///
  /// A zero-mutation branch has a finite, maximal density at the boundary, so the neg-log is a
  /// straight line. The law is that exact line through the boundary point `(t_hard, y_hard)` and the
  /// grid edge -- two points, no regression: `slope = (y_edge - y_hard) / (t_edge - t_hard)`.
  /// `y_hard = -ln p(t_hard)` is the boundary neg-log ordinate, supplied by the producer (the grid
  /// samples strictly inside `min_bl > 0`, so no grid point sits at `t_hard`).
  ///
  /// `Side` selects the grid end nearest the boundary. Returns an error on an empty grid, a
  /// non-finite edge ordinate, or a non-finite slope.
  pub fn fit_linear(grid_fn: &GridFn<f64>, t_hard: f64, side: Side, y_hard: f64) -> Result<Self, Report> {
    let (t_edge, y_edge) = edge_ordinate(grid_fn, side)?;
    let slope = (y_edge - y_hard) / (t_edge - t_hard);
    if !slope.is_finite() {
      return make_error!(
        "Hard-boundary linear fit on the {side:?} side produced a non-finite slope from grid edge \
         (t={t_edge}, y={y_edge}) and boundary (t={t_hard}, y={y_hard})"
      );
    }
    Ok(HardApproachLaw {
      t_hard,
      shape: Approach::Finite { slope },
    })
  }

  /// Fit the power-law approach for a divergent boundary density (`n >= 1`, or an indel) over
  /// `[t_hard, t_edge)`.
  ///
  /// The density vanishes at the boundary, so `-ln p` diverges. Fit `b >= 0` from the innermost
  /// `n_fit` points by one log-distance regression on `(ln|t - t_hard|, y)`: `y = a - b*ln|dt|` is
  /// linear in `ln|dt|` with slope `-b`. The intercept `a` is discarded (edge-relative). A
  /// wrong-sign fit clamps to `b = 0`, which yields a flat [`Finite`](Approach::Finite) shape rather
  /// than a divergent one that would fabricate a boundary singularity the data does not support.
  ///
  /// `Side` selects the grid end nearest the boundary. Returns an error on an empty grid, a
  /// non-finite edge ordinate, fewer than two finite innermost points, or a non-finite result.
  pub fn fit_log_power_law(grid_fn: &GridFn<f64>, t_hard: f64, side: Side, n_fit: usize) -> Result<Self, Report> {
    let n = grid_fn.n_points();
    // Validate that the grid can anchor an edge-relative law at evaluation time; the regression
    // below reads interior points, but a non-finite or absent edge cannot anchor the live law.
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
    // A wrong-sign regression clamps to `b = 0`, which is the flat finite shape, not a divergence.
    let shape = if b > 0.0 {
      Approach::Divergent { b }
    } else {
      Approach::Finite { slope: 0.0 }
    };
    Ok(HardApproachLaw { t_hard, shape })
  }

  /// Evaluate the approach law in neg-log at `t`, anchored on the live grid `edge`.
  ///
  /// `edge` is the grid's current edge (its coordinate `edge.t` and stored neg-log ordinate
  /// `edge.y`). Each regime is its own closed form, so the law meets the grid continuously at the
  /// edge:
  ///
  /// - [`Finite`](Approach::Finite): the line `edge.y + slope*(t - edge.t)`, finite everywhere
  ///   including at `t == t_hard`, where it carries the boundary mode.
  /// - [`Divergent`](Approach::Divergent): the power-law
  ///   `edge.y - b*ln(|t - t_hard| / |edge.t - t_hard|)`, returning `+inf` at `t == t_hard` because
  ///   the density is zero at the boundary.
  pub fn eval(&self, edge: GridEdge, t: f64) -> f64 {
    match self.shape {
      Approach::Finite { slope } => edge.y + slope * (t - edge.t),
      Approach::Divergent { b } => {
        let dt = (t - self.t_hard).abs();
        if dt == 0.0 {
          return f64::INFINITY;
        }
        let dt_edge = (edge.t - self.t_hard).abs();
        edge.y - b * (dt / dt_edge).ln()
      },
    }
  }

  /// Probability mass in the approach region between `t_hard` and the grid edge, in plain
  /// probability.
  ///
  /// Closed form of `integral_{t_hard}^{edge.t} exp(-y(t)) dt` with the edge probability recovered
  /// from the anchor's stored neg-log ordinate as `p_edge = exp(-edge.y)`:
  ///
  /// - [`Finite`](Approach::Finite), `slope != 0`: `p_edge * (exp(slope*dt_edge) - 1) / slope`.
  /// - [`Finite`](Approach::Finite), `slope == 0`: `p_edge * dt_edge`.
  /// - [`Divergent`](Approach::Divergent): `p_edge * dt_edge / (b + 1)`.
  pub fn mass(&self, edge: GridEdge) -> f64 {
    let dt_edge = (edge.t - self.t_hard).abs();
    let p_edge = (-edge.y).exp();

    match self.shape {
      Approach::Finite { slope } if slope != 0.0 => p_edge * (slope * dt_edge).exp_m1() / slope,
      Approach::Finite { .. } => p_edge * dt_edge,
      Approach::Divergent { b } => p_edge * dt_edge / (b + 1.0),
    }
  }

  /// Whether the density is finite at the boundary (the [`Finite`](Approach::Finite) regime).
  ///
  /// `true` means the neg-log reaches a finite value at `t_hard` and the mode sits on the boundary;
  /// `false` means the density vanishes there ([`Divergent`](Approach::Divergent), `y -> +inf`).
  /// Consumers use this to ask the semantic question instead of inspecting the raw exponent.
  #[must_use]
  pub fn is_finite_at_boundary(&self) -> bool {
    matches!(self.shape, Approach::Finite { .. })
  }

  /// Scale the approach law when every neg-log ordinate is multiplied by `factor` (`p -> p^factor`).
  ///
  /// The regime's shape parameter scales by `factor`; the regime is preserved for `factor > 0`.
  #[must_use]
  pub fn scale(&self, factor: f64) -> HardApproachLaw {
    let shape = match self.shape {
      Approach::Finite { slope } => Approach::Finite { slope: slope * factor },
      Approach::Divergent { b } => Approach::Divergent { b: b * factor },
    };
    HardApproachLaw {
      t_hard: self.t_hard,
      shape,
    }
  }

  /// Transform the approach law when the argument is negated: `f(x) -> f(-x)`.
  ///
  /// The boundary moves to `-t_hard`; the divergent exponent `b` is unchanged and the finite slope
  /// flips sign.
  #[must_use]
  pub fn negate_arg(&self) -> Self {
    let shape = match self.shape {
      Approach::Finite { slope } => Approach::Finite { slope: -slope },
      Approach::Divergent { b } => Approach::Divergent { b },
    };
    HardApproachLaw {
      t_hard: -self.t_hard,
      shape,
    }
  }
}

/// The neg-log approach shape near a hard boundary, in one of two mutually exclusive regimes.
///
/// The regime is decided once, by the producer that owns the likelihood, and stored in the type so
/// no consumer re-derives it from float values. The branch-length likelihood `p(t) ~ t^n *
/// exp(-mu*t)` (`n` mutations) gives both regimes directly.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Approach {
  /// Finite boundary density (`n = 0`): the neg-log is the straight line
  /// `y_edge + slope*(t - t_edge)`, finite and maximal at the boundary, which the linear term
  /// carries. `slope == 0` is the flat sub-case.
  Finite { slope: f64 },
  /// Divergent boundary density (`n >= 1`, or an indel): the density vanishes at the boundary, so
  /// the neg-log is the power law `y_edge - b*ln|dt/dt_edge|` with `b > 0` and `y -> +inf` at
  /// `t_hard`. `b` is a fitted continuous approximation of the mutation count, not an integer.
  Divergent { b: f64 },
}

/// Which side of the grid a hard boundary is on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
  Left,
  Right,
}

/// Grid edge nearest the boundary as `(t_edge, y_edge)`.
///
/// Returns an error on an empty grid or a non-finite edge ordinate, neither of which can anchor an
/// edge-relative law.
fn edge_ordinate(grid_fn: &GridFn<f64>, side: Side) -> Result<(f64, f64), Report> {
  let n = grid_fn.n_points();
  if n == 0 {
    return make_error!("Cannot anchor an edge-relative law on an empty grid");
  }
  let edge_idx = match side {
    Side::Left => 0,
    Side::Right => n - 1,
  };
  let t_edge = grid_fn.grid().x_at(edge_idx);
  let y_edge = grid_fn.y()[edge_idx];
  if !y_edge.is_finite() {
    return make_error!(
      "Grid {side:?} edge ordinate is non-finite (y={y_edge} at t={t_edge}); cannot anchor an \
       edge-relative law"
    );
  }
  Ok((t_edge, y_edge))
}
