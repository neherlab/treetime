use crate::GridFn;
use eyre::Report;
use serde::{Deserialize, Serialize};
use treetime_utils::make_error;

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
/// right) and `t_edge` its coordinate. The law is *edge-relative*: it stores only the shape and
/// reads the live grid edge on evaluation, exactly like [`SoftTailLaw`](crate::SoftTailLaw).
/// Reading the live edge keeps the law valid across re-windowing and resampling with no refit and
/// no absolute anchor to keep in sync: a vertical shift of every ordinate (peak-normalization) is
/// absorbed through `y_edge`, and a scale (`p -> p^factor`) scales the shape.
///
/// The boundary location `t_hard` is an immovable physical fact (e.g. `t = 0` for branch lengths)
/// and is stored as an absolute coordinate. The shape is one of three [regimes](Approach), decided once at
/// construction and never re-derived by a consumer: every evaluator matches [`Approach`] instead
/// of inspecting float values, so `eval` and `mass` cannot disagree about which regime a law is.
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
    Ok(HardApproachLaw::from_terms(t_hard, 0.0, slope))
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

    let neg_b_raw = least_squares_slope(&xs, &ys);
    let b = (-neg_b_raw).max(0.0);
    if !b.is_finite() {
      return make_error!("Hard-boundary power-law fit on the {side:?} side produced a non-finite exponent");
    }
    Ok(HardApproachLaw::from_terms(t_hard, b, 0.0))
  }

  /// Build a law from raw `b` and `slope` terms, classifying the [regime](Approach) once.
  ///
  /// This is the single point where the two shape parameters are mapped to a regime: `b > 0` with a
  /// non-zero slope is [`Combined`](Approach::Combined), `b > 0` alone is
  /// [`Divergent`](Approach::Divergent), and everything else (including a clamped `b = 0`) is
  /// [`Finite`](Approach::Finite). Producers (`fit_*`, `compose_multiply`, `scale`, `negate_arg`)
  /// go through here; consumers (`eval`, `mass`) match the resulting regime and never re-derive it.
  pub(crate) fn from_terms(t_hard: f64, b: f64, slope: f64) -> Self {
    let shape = if b > 0.0 {
      if slope != 0.0 {
        Approach::Combined { b, slope }
      } else {
        Approach::Divergent { b }
      }
    } else {
      Approach::Finite { slope }
    };
    HardApproachLaw { t_hard, shape }
  }

  /// The shape terms of this law, with the regime's implied zero filled in.
  fn terms(&self) -> ShapeTerms {
    match self.shape {
      Approach::Finite { slope } => ShapeTerms { b: 0.0, slope },
      Approach::Divergent { b } => ShapeTerms { b, slope: 0.0 },
      Approach::Combined { b, slope } => ShapeTerms { b, slope },
    }
  }

  /// Evaluate the approach law in neg-log at `t`, anchored on the live grid edge.
  ///
  /// `y_edge` is the grid's stored neg-log edge ordinate and `t_edge` its coordinate. Each regime is
  /// its own closed form, so the law meets the grid continuously at the edge:
  ///
  /// - [`Finite`](Approach::Finite): the line `y_edge + slope*(t - t_edge)`, finite everywhere
  ///   including at `t == t_hard`, where it carries the boundary mode.
  /// - [`Divergent`](Approach::Divergent) and [`Combined`](Approach::Combined): the power-law
  ///   `y_edge - b*ln(|t - t_hard| / |t_edge - t_hard|)` (plus the `slope` term for `Combined`),
  ///   returning `+inf` at `t == t_hard` because the density is zero at the boundary.
  pub fn eval(&self, y_edge: f64, t_edge: f64, t: f64) -> f64 {
    match self.shape {
      Approach::Finite { slope } => y_edge + slope * (t - t_edge),
      Approach::Divergent { b } => {
        let dt = (t - self.t_hard).abs();
        if dt == 0.0 {
          return f64::INFINITY;
        }
        let dt_edge = (t_edge - self.t_hard).abs();
        y_edge - b * (dt / dt_edge).ln()
      },
      Approach::Combined { b, slope } => {
        let dt = (t - self.t_hard).abs();
        if dt == 0.0 {
          return f64::INFINITY;
        }
        let dt_edge = (t_edge - self.t_hard).abs();
        y_edge - b * (dt / dt_edge).ln() + slope * (t - t_edge)
      },
    }
  }

  /// Probability mass in the approach region between `t_hard` and the grid edge, in plain
  /// probability.
  ///
  /// Closed form of `integral_{t_hard}^{t_edge} exp(-y(t)) dt` with the edge probability recovered
  /// from its stored neg-log ordinate as `p_edge = exp(-y_edge)`:
  ///
  /// - [`Finite`](Approach::Finite), `slope != 0`: `p_edge * (exp(slope*dt_edge) - 1) / slope`.
  /// - [`Finite`](Approach::Finite), `slope == 0`: `p_edge * dt_edge`.
  /// - [`Divergent`](Approach::Divergent): `p_edge * dt_edge / (b + 1)`.
  /// - [`Combined`](Approach::Combined): `p_edge * dt_edge / (b + 1)` -- the `slope` correction is
  ///   dropped, so this mass does *not* match the integral of [`eval`](Self::eval) for a `Combined`
  ///   law. See `kb/issues/M-grid-hard-approach-combined-mass-drops-slope.md` for the exact
  ///   incomplete-gamma form pending a decision.
  pub fn mass(&self, y_edge: f64, t_edge: f64) -> f64 {
    let dt_edge = (t_edge - self.t_hard).abs();
    let p_edge = (-y_edge).exp();

    match self.shape {
      Approach::Finite { slope } if slope != 0.0 => p_edge * (slope * dt_edge).exp_m1() / slope,
      Approach::Finite { .. } => p_edge * dt_edge,
      Approach::Divergent { b } | Approach::Combined { b, .. } => p_edge * dt_edge / (b + 1.0),
    }
  }

  /// Whether the density is finite at the boundary (the [`Finite`](Approach::Finite) regime).
  ///
  /// `true` means the neg-log reaches a finite value at `t_hard` and the mode sits on the boundary;
  /// `false` means the density vanishes there (`Divergent`/`Combined`, `y -> +inf`). Consumers use
  /// this to ask the semantic question instead of inspecting the raw exponent.
  #[must_use]
  pub fn is_finite_at_boundary(&self) -> bool {
    matches!(self.shape, Approach::Finite { .. })
  }

  /// Compose two approach laws under multiplication: the shape terms add.
  ///
  /// Multiplication is addition in neg-log space, so the product law's exponent and slope are the
  /// sums of the operand terms, re-classified into a regime. A product of a divergent and a finite
  /// message therefore carries both a `b > 0` and a `slope != 0` term and lands in
  /// [`Combined`](Approach::Combined). Both laws must share the same `t_hard`, and both are evaluated
  /// against the same live grid edge, so there is no anchor to reconcile.
  #[must_use]
  pub fn compose_multiply(&self, other: &HardApproachLaw) -> HardApproachLaw {
    let a = self.terms();
    let b = other.terms();
    HardApproachLaw::from_terms(self.t_hard, a.b + b.b, a.slope + b.slope)
  }

  /// Scale the approach law when every neg-log ordinate is multiplied by `factor` (`p -> p^factor`).
  ///
  /// Both shape terms scale by `factor`; the regime is preserved for `factor > 0`.
  #[must_use]
  pub fn scale(&self, factor: f64) -> HardApproachLaw {
    let t = self.terms();
    HardApproachLaw::from_terms(self.t_hard, t.b * factor, t.slope * factor)
  }

  /// Transform the approach law when the argument is negated: `f(x) -> f(-x)`.
  ///
  /// The boundary moves to `-t_hard`; the exponent `b` is unchanged and the linear slope flips sign.
  #[must_use]
  pub fn negate_arg(&self) -> Self {
    let t = self.terms();
    HardApproachLaw::from_terms(-self.t_hard, t.b, -t.slope)
  }
}

/// The two raw shape terms of a hard-approach law before classification into an [`Approach`] regime.
struct ShapeTerms {
  b: f64,
  slope: f64,
}

/// The neg-log approach shape near a hard boundary, in one of three mutually exclusive regimes.
///
/// The regime is decided once, by the producer that owns the likelihood, and stored in the type so
/// no consumer re-derives it from float values. The branch-length likelihood `p(t) ~ t^n *
/// exp(-mu*t)` (`n` mutations) gives the first two regimes directly; multiplication of messages
/// (addition in neg-log space) can combine them into the third.
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
  /// Product of a divergent and a finite message: both terms are present (`b > 0`, `slope != 0`).
  /// This regime exists because composition adds shape terms and is not closed over the first two
  /// alone.
  Combined { b: f64, slope: f64 },
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

/// Slope of the simple least-squares line `y = slope * x + intercept`.
///
/// The intercept is not returned: both call sites (the power-law exponent fit and the soft-tail
/// slope fit) are edge-relative and discard it.
pub(crate) fn least_squares_slope(xs: &[f64], ys: &[f64]) -> f64 {
  let n = xs.len() as f64;
  let sum_x: f64 = xs.iter().sum();
  let sum_y: f64 = ys.iter().sum();
  let sum_xx: f64 = xs.iter().map(|x| x * x).sum();
  let sum_xy: f64 = xs.iter().zip(ys).map(|(x, y)| x * y).sum();

  let sum_x_sq = sum_x * sum_x;
  let denom = n * sum_xx - sum_x_sq;
  if denom.abs() < 1e-30 {
    return 0.0;
  }

  (n * sum_xy - sum_x * sum_y) / denom
}
