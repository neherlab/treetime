use crate::hard_approach_law::Side;
use crate::{GridEdge, GridFn};
use eyre::Report;
use serde::{Deserialize, Serialize};
use treetime_utils::least_squares::LineFit;
use treetime_utils::make_error;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SoftTailLaw {
  pub slope: f64,
}

impl SoftTailLaw {
  pub fn fit(grid_fn: &GridFn<f64>, side: Side, n_fit: usize) -> Result<Self, Report> {
    let n = grid_fn.n_points();
    let n_fit = n_fit.min(n);

    let t_edge = match side {
      Side::Left => grid_fn.grid().x_at(0),
      Side::Right => grid_fn.grid().x_at(n - 1),
    };

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

    let slope = match side {
      Side::Left => slope_raw.min(0.0),
      Side::Right => slope_raw.max(0.0),
    };

    Ok(SoftTailLaw { slope })
  }

  pub fn eval(&self, edge: GridEdge, t: f64) -> f64 {
    edge.y + self.slope * (t - edge.t)
  }

  pub fn mass(&self, y_edge: f64) -> f64 {
    (-y_edge).exp() / self.slope.abs()
  }

  #[must_use]
  pub fn compose_multiply(&self, other: &SoftTailLaw) -> SoftTailLaw {
    SoftTailLaw {
      slope: self.slope + other.slope,
    }
  }

  #[must_use]
  pub fn negate_arg(&self) -> Self {
    SoftTailLaw { slope: -self.slope }
  }
}
