use crate::{GridEdge, GridFn};
use eyre::Report;
use serde::{Deserialize, Serialize};
use treetime_utils::least_squares::LineFit;
use treetime_utils::make_error;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct HardApproachLaw {
  pub t_hard: f64,
  pub b: f64,
}

impl HardApproachLaw {
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

  pub(crate) fn eval(&self, edge: GridEdge, t: f64) -> f64 {
    let dt = (t - self.t_hard).abs();
    if dt == 0.0 {
      return if self.b == 0.0 { edge.y } else { f64::INFINITY };
    }
    let dt_edge = (edge.t - self.t_hard).abs();
    edge.y - self.b * (dt / dt_edge).ln()
  }

  pub fn mass(&self, edge: GridEdge) -> f64 {
    let dt_edge = (edge.t - self.t_hard).abs();
    let p_edge = (-edge.y).exp();
    p_edge * dt_edge / (self.b + 1.0)
  }

  #[must_use]
  pub fn negate_arg(&self) -> Self {
    HardApproachLaw {
      t_hard: -self.t_hard,
      b: self.b,
    }
  }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Side {
  Left,
  Right,
}
