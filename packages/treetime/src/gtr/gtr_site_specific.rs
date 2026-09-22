use crate::gtr::gtr::eig_single_site;
use crate::{make_error, make_report};
use eyre::Report;
use ndarray::prelude::*;
use ndarray::{Array3, Array4};
use rand::Rng;
use rand_distr::Gamma;
use treetime_utils::array::ndarray::clamp_min;

#[derive(Clone, Debug)]
pub struct GTRSiteSpecificParams {
  pub n_states: usize,
  pub seq_len: usize,
  pub mu: Array1<f64>,
  pub W: Option<Array2<f64>>,
  pub pi: Array2<f64>,
  pub approximate: bool,
}

#[derive(Clone, Debug)]
#[allow(clippy::partial_pub_fields)]
pub struct GTRSiteSpecific {
  pub seq_len: usize,
  pub mu: Array1<f64>,
  pub W: Array2<f64>,
  pub pi: Array2<f64>,
  pub eigvals: Array2<f64>,
  v: Array3<f64>,
  v_inv: Array3<f64>,
  interpolator: Option<ExpQtInterpolator>,
}

impl GTRSiteSpecific {
  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
  pub(crate) fn new(
    GTRSiteSpecificParams {
      n_states,
      seq_len,
      mu,
      W,
      pi,
      approximate,
    }: GTRSiteSpecificParams,
  ) -> Result<Self, Report> {
    if mu.len() != seq_len {
      return make_error!("mu length {} does not match seq_len {seq_len}", mu.len());
    }
    if pi.dim() != (n_states, seq_len) {
      return make_error!(
        "pi shape {:?} does not match expected ({n_states}, {seq_len})",
        pi.dim()
      );
    }
    if let Some(W) = &W {
      if W.dim() != (n_states, n_states) {
        return make_error!("W shape {:?} does not match expected ({n_states}, {n_states})", W.dim());
      }
    }
    for a in 0..seq_len {
      if mu[a] < 0.0 {
        return make_error!("Site {a} has negative substitution rate mu: {}", mu[a]);
      }
    }

    let W = {
      let W = W.unwrap_or_else(|| {
        let mut W = Array2::<f64>::ones([n_states, n_states]);
        W.diag_mut().fill(0.0);
        W
      });
      let mut W = 0.5 * (&W.view() + &W.t());
      W.diag_mut().fill(0.0);
      W
    };

    let pi = {
      let col_sums = pi.sum_axis(Axis(0));
      for a in 0..seq_len {
        if col_sums[a] <= 0.0 {
          return make_error!("Site {a} has non-positive pi column sum: {}", col_sums[a]);
        }
      }
      let pi = &pi / &col_sums;
      for a in 0..seq_len {
        for i in 0..n_states {
          if pi[[i, a]] <= 0.0 {
            return make_error!(
              "Site {a}, state {i} has non-positive pi after normalization: {}",
              pi[[i, a]]
            );
          }
        }
      }
      pi
    };

    let mut mu = mu;
    let W = {
      let mut total_avg = 0.0;
      for a in 0..seq_len {
        let pi_a = pi.column(a);
        let rate_a = pi_a.dot(&W.dot(&pi_a));
        total_avg += rate_a;
      }
      total_avg /= seq_len as f64;
      if total_avg <= 0.0 {
        return make_error!("Average substitution rate is non-positive: {total_avg}");
      }
      mu *= total_avg;
      W / total_avg
    };

    let mut eigvals = Array2::zeros((n_states, seq_len));
    let mut v = Array3::zeros((n_states, n_states, seq_len));
    let mut v_inv = Array3::zeros((n_states, n_states, seq_len));

    for a in 0..seq_len {
      let (ev, evec, evec_inv) = eig_single_site(&W, pi.column(a))?;
      eigvals.column_mut(a).assign(&ev);
      v.slice_mut(s![.., .., a]).assign(&evec);
      v_inv.slice_mut(s![.., .., a]).assign(&evec_inv);
    }

    let mut model = Self {
      seq_len,
      mu,
      W,
      pi,
      eigvals,
      v,
      v_inv,
      interpolator: None,
    };

    if approximate {
      model.build_interpolator();
    }

    Ok(model)
  }

  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
  pub fn random(
    n_states: usize,
    seq_len: usize,
    avg_mu: f64,
    pi_dirichlet_alpha: f64,
    W_dirichlet_alpha: f64,
    mu_gamma_alpha: f64,
    rng: &mut impl Rng,
  ) -> Result<Self, Report> {
    let pi = {
      let mut pi = Array2::zeros((n_states, seq_len));
      if pi_dirichlet_alpha > 0.0 {
        let gamma = Gamma::new(pi_dirichlet_alpha, 1.0).map_err(|e| make_report!("{e}"))?;
        for a in 0..seq_len {
          for i in 0..n_states {
            pi[[i, a]] = rng.sample(gamma);
          }
        }
      } else {
        pi.fill(1.0);
      }
      let col_sums = pi.sum_axis(Axis(0));
      &pi / &col_sums
    };

    let W = {
      let mut W = Array2::zeros((n_states, n_states));
      if W_dirichlet_alpha > 0.0 {
        let gamma = Gamma::new(W_dirichlet_alpha, 1.0).map_err(|e| make_report!("{e}"))?;
        for i in 0..n_states {
          for j in 0..i {
            let val: f64 = rng.sample(gamma);
            W[[i, j]] = val;
            W[[j, i]] = val;
          }
        }
      } else {
        W.fill(1.0);
        W.diag_mut().fill(0.0);
      }
      W
    };

    let mu = {
      let mut mu = Array1::zeros(seq_len);
      if mu_gamma_alpha > 0.0 {
        let gamma = Gamma::new(mu_gamma_alpha, 1.0).map_err(|e| make_report!("{e}"))?;
        for a in 0..seq_len {
          mu[a] = rng.sample(gamma);
        }
      } else {
        mu.fill(1.0);
      }
      mu
    };

    let mut model = Self::new(GTRSiteSpecificParams {
      n_states,
      seq_len,
      mu,
      W: Some(W),
      pi,
      approximate: false,
    })?;

    let mean_rate = model.average_rate().sum() / seq_len as f64;
    if mean_rate > 0.0 {
      model.mu *= avg_mu / mean_rate;
    }

    Ok(model)
  }

  pub(crate) fn average_rate(&self) -> Array1<f64> {
    let mut rates = Array1::zeros(self.seq_len);
    for a in 0..self.seq_len {
      let pi_a = self.pi.column(a);
      rates[a] = self.mu[a] * pi_a.dot(&self.W.dot(&pi_a));
    }
    rates
  }

  pub fn expQt(&self, t: f64) -> Result<Array3<f64>, Report> {
    if t < 0.0 {
      return make_error!("Branch length t must be non-negative, got {t}");
    }
    if let Some(interp) = &self.interpolator {
      if t * interp.rate_scale < ExpQtInterpolator::MAX_INTERP_RANGE {
        return Ok(interp.interpolate(t));
      }
    }
    Ok(self.expQt_raw(t))
  }

  pub(crate) fn expQt_raw(&self, t: f64) -> Array3<f64> {
    let n = self.eigvals.nrows();
    let mut result = Array3::zeros((n, n, self.seq_len));

    for a in 0..self.seq_len {
      let e_lambda_t: Array1<f64> = (&self.eigvals.column(a) * self.mu[a] * t).mapv(f64::exp);
      let v_a = self.v.slice(s![.., .., a]);
      let v_inv_a = self.v_inv.slice(s![.., .., a]);
      let scaled_v = &v_a * &e_lambda_t;
      let p_a = scaled_v.dot(&v_inv_a);
      result.slice_mut(s![.., .., a]).assign(&clamp_min(&p_a, 0.0));
    }

    result
  }

  pub fn propagate_profile(&self, profile: &Array2<f64>, t: f64, return_log: bool) -> Result<Array2<f64>, Report> {
    let qt = self.expQt(t)?;
    let mut result = Array2::zeros(profile.dim());

    for a in 0..self.seq_len {
      let qt_a = qt.slice(s![.., .., a]);
      result.row_mut(a).assign(&profile.row(a).dot(&qt_a));
    }

    if return_log {
      result.mapv_inplace(f64::ln);
    }
    Ok(result)
  }

  pub fn evolve(&self, profile: &Array2<f64>, t: f64, return_log: bool) -> Result<Array2<f64>, Report> {
    let qt = self.expQt(t)?;
    let mut result = Array2::zeros(profile.dim());

    for a in 0..self.seq_len {
      let qt_a = qt.slice(s![.., .., a]);
      result.row_mut(a).assign(&profile.row(a).dot(&qt_a.t()));
    }

    if return_log {
      result.mapv_inplace(f64::ln);
    }
    Ok(result)
  }

  pub fn Q(&self) -> Array3<f64> {
    let n = self.pi.nrows();
    let mut result = Array3::zeros((n, n, self.seq_len));
    for a in 0..self.seq_len {
      let pi_a = self.pi.column(a);
      let mut q_a = (&self.W * &pi_a).t().to_owned();
      let diag = -q_a.sum_axis(Axis(0));
      q_a.diag_mut().assign(&diag);
      result.slice_mut(s![.., .., a]).assign(&q_a);
    }
    result
  }

  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
  fn build_interpolator(&mut self) {
    let avg_rates = self.average_rate();
    let rate_scale = (avg_rates.sum() / self.seq_len as f64).max(1e-10);
    let inv_rate = 1.0 / rate_scale;

    let mut t_grid = Vec::new();
    for &v in &linspace(0.0, 0.1, 11)[..10] {
      t_grid.push(v * inv_rate);
    }
    for &v in &linspace(0.1, 1.0, 21)[..20] {
      t_grid.push(v * inv_rate);
    }
    for &v in &linspace(1.0, 5.0, 21)[..20] {
      t_grid.push(v * inv_rate);
    }
    for &v in &linspace(5.0, 10.0, 11) {
      t_grid.push(v * inv_rate);
    }
    let t_grid = Array1::from_vec(t_grid);

    let n_t = t_grid.len();
    let n = self.eigvals.nrows();
    let mut data = Array4::zeros((n_t, n, n, self.seq_len));
    for (idx, &t) in t_grid.iter().enumerate() {
      let qt = self.expQt_raw(t);
      data.slice_mut(s![idx, .., .., ..]).assign(&qt);
    }

    self.interpolator = Some(ExpQtInterpolator {
      t_grid,
      data,
      rate_scale,
    });
  }
}

#[derive(Clone, Debug)]
pub struct ExpQtInterpolator {
  pub t_grid: Array1<f64>,
  pub data: Array4<f64>,
  pub rate_scale: f64,
}

impl ExpQtInterpolator {
  pub const MAX_INTERP_RANGE: f64 = 10.0;

  pub fn interpolate(&self, t: f64) -> Array3<f64> {
    let n = self.t_grid.len();
    debug_assert!(n >= 2);

    let t = t.clamp(self.t_grid[0], self.t_grid[n - 1]);

    let left = self
      .t_grid
      .iter()
      .position(|&v| v > t)
      .map_or(n - 2, |idx| idx.saturating_sub(1).min(n - 2));

    let t_left = self.t_grid[left];
    let t_right = self.t_grid[left + 1];
    let alpha = if (t_right - t_left).abs() < 1e-15 {
      0.0
    } else {
      (t - t_left) / (t_right - t_left)
    };

    let left_slice = self.data.slice(s![left, .., .., ..]);
    let right_slice = self.data.slice(s![left + 1, .., .., ..]);

    &left_slice * (1.0 - alpha) + &right_slice * alpha
  }
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
fn linspace(start: f64, end: f64, n: usize) -> Vec<f64> {
  if n <= 1 {
    return vec![start];
  }
  let step = (end - start) / (n - 1) as f64;
  (0..n).map(|i| start + i as f64 * step).collect()
}
