use crate::gtr::gtr::eig_single_site;
use crate::make_error;
use eyre::Report;
use ndarray::prelude::*;
use ndarray::{Array3, Array4};
use rand::Rng;
use rand_distr::Gamma;
use treetime_utils::array::ndarray::clamp_min;

/// Parameters for constructing a site-specific GTR model.
#[derive(Clone, Debug)]
pub struct GTRSiteSpecificParams {
  /// Number of character states (e.g. 4 for nucleotides, 20 for amino acids).
  pub n_states: usize,
  /// Sequence length (number of alignment positions).
  pub seq_len: usize,
  /// Per-site substitution rate. Shape: [seq_len].
  pub mu: Array1<f64>,
  /// Symmetric exchangeability matrix (shared across all sites). Shape: [n_states, n_states].
  /// If None, defaults to all-ones (equal rates). Will be symmetrized and diagonal zeroed.
  pub W: Option<Array2<f64>>,
  /// Per-site equilibrium frequencies. Shape: [n_states, seq_len].
  /// Each column sums to 1. This is what makes the model "site-specific":
  /// different sites have different stationary distributions, requiring
  /// per-site eigendecomposition of the rate matrix.
  pub pi: Array2<f64>,
  /// Use linear interpolation of pre-computed exp(Qt) for fast evaluation.
  pub approximate: bool,
}

/// Site-specific General Time-Reversible model of character evolution.
///
/// Extension of GTR where equilibrium frequencies vary per alignment site. Because the
/// rate matrix Q depends on pi, and pi differs across sites, each site requires its own
/// eigendecomposition. This makes the transition probability matrix P(t) = exp(Q*mu*t)
/// site-specific: a 3D array [n_states, n_states, seq_len] rather than a single 2D matrix.
///
/// The shared exchangeability matrix W encodes relative rates between states.
/// Per-site pi values then produce per-site rate matrices Q_a = f(W, pi_a).
///
/// For computational efficiency, an optional interpolation mode pre-computes exp(Q*t)
/// on a grid of t values and uses linear interpolation during tree traversal.
///
/// Reference: Puller et al (2020). "TreeTime: Maximum-likelihood phylodynamic analysis."
/// Virus Evolution, 6(1), veaa066.
#[derive(Clone, Debug)]
#[allow(clippy::partial_pub_fields)]
pub struct GTRSiteSpecific {
  /// Sequence length (number of alignment sites).
  pub seq_len: usize,
  /// Per-site substitution rates. Shape: [seq_len].
  pub mu: Array1<f64>,
  /// Symmetric exchangeability matrix (shared, zero diagonal, normalized). Shape: [n_states, n_states].
  pub W: Array2<f64>,
  /// Per-site equilibrium frequencies. Shape: [n_states, seq_len]. Each column sums to 1.
  pub pi: Array2<f64>,
  /// Per-site eigenvalues. Shape: [n_states, seq_len].
  /// All values are non-positive; exactly one zero eigenvalue per site.
  pub eigvals: Array2<f64>,
  /// Per-site right transformation matrices. Shape: [n_states, n_states, seq_len].
  pub v: Array3<f64>,
  /// Per-site left transformation matrices. Shape: [n_states, n_states, seq_len].
  pub v_inv: Array3<f64>,
  /// Pre-computed interpolation data for fast exp(Qt) evaluation.
  interpolator: Option<ExpQtInterpolator>,
}

impl GTRSiteSpecific {
  /// Construct a new site-specific GTR model.
  ///
  /// Performs per-site eigendecomposition of the rate matrix Q_a = f(W, pi_a) for each
  /// alignment position a. Optionally builds an interpolation table for fast exp(Qt).
  ///
  /// # Processing steps
  ///
  /// 1. Symmetrize W and zero diagonal (shared across sites).
  /// 2. Normalize pi columns to sum to 1.
  /// 3. Normalize W by average rate across all sites, absorb scaling into mu.
  /// 4. Eigendecompose Q_a for each site a, producing per-site eigvals, v, v_inv.
  /// 5. Optionally pre-compute exp(Qt) interpolation table.
  pub fn new(
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

    // Symmetrize W and zero diagonal
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

    // Validate and normalize pi columns to sum to 1.
    // The similarity transform in eig_single_site divides by sqrt(pi), so all
    // entries must be strictly positive. A zero column sum or zero entry would
    // produce NaN/Inf in the eigendecomposition.
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

    // Normalize W by average rate across all sites.
    // average_rate_a = pi_a^T * W * pi_a for each site a.
    // Global average = mean over sites.
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

    // Per-site eigendecomposition
    let mut eigvals = Array2::zeros((n_states, seq_len));
    let mut v = Array3::zeros((n_states, n_states, seq_len));
    let mut v_inv = Array3::zeros((n_states, n_states, seq_len));

    for a in 0..seq_len {
      let pi_a = pi.column(a).to_owned();
      let (ev, evec, evec_inv) = eig_single_site(&W, &pi_a)?;
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

  /// Generate a random site-specific GTR model from prior distributions.
  ///
  /// Samples per-site equilibrium frequencies from Dirichlet (via Gamma),
  /// shared exchangeability matrix W from Gamma, and per-site rates from Gamma.
  /// Matches v0's `GTR_site_specific.random()`.
  ///
  /// # Arguments
  ///
  /// * `n_states` - Number of character states (e.g. 4 for nucleotides).
  /// * `seq_len` - Number of alignment sites.
  /// * `avg_mu` - Target average substitution rate.
  /// * `pi_dirichlet_alpha` - Dirichlet concentration for per-site equilibrium frequencies.
  ///   0 produces uniform pi.
  /// * `W_dirichlet_alpha` - Gamma shape for exchangeability matrix entries. 0 produces uniform W.
  /// * `mu_gamma_alpha` - Gamma shape for per-site rates. 0 produces uniform mu.
  /// * `rng` - Random number generator.
  pub fn random(
    n_states: usize,
    seq_len: usize,
    avg_mu: f64,
    pi_dirichlet_alpha: f64,
    W_dirichlet_alpha: f64,
    mu_gamma_alpha: f64,
    rng: &mut impl Rng,
  ) -> Result<Self, Report> {
    // Dirichlet-distributed pi: sample Gamma per entry, then L1-normalize columns
    let pi = {
      let mut pi = Array2::zeros((n_states, seq_len));
      if pi_dirichlet_alpha > 0.0 {
        let gamma = Gamma::new(pi_dirichlet_alpha, 1.0).map_err(|e| eyre::eyre!("{e}"))?;
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

    // Symmetric exchangeability matrix from Gamma-distributed lower triangle
    let W = {
      let mut W = Array2::zeros((n_states, n_states));
      if W_dirichlet_alpha > 0.0 {
        let gamma = Gamma::new(W_dirichlet_alpha, 1.0).map_err(|e| eyre::eyre!("{e}"))?;
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

    // Per-site rates from Gamma distribution
    let mu = {
      let mut mu = Array1::zeros(seq_len);
      if mu_gamma_alpha > 0.0 {
        let gamma = Gamma::new(mu_gamma_alpha, 1.0).map_err(|e| eyre::eyre!("{e}"))?;
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

    // Scale mu so that mean average_rate equals avg_mu
    let mean_rate = model.average_rate().sum() / seq_len as f64;
    if mean_rate > 0.0 {
      model.mu *= avg_mu / mean_rate;
    }

    Ok(model)
  }

  /// Compute the per-site average substitution rate.
  ///
  /// Returns mu_a * pi_a^T * W * pi_a for each site a.
  /// Shape: [seq_len].
  pub fn average_rate(&self) -> Array1<f64> {
    let mut rates = Array1::zeros(self.seq_len);
    for a in 0..self.seq_len {
      let pi_a = self.pi.column(a);
      rates[a] = self.mu[a] * pi_a.dot(&self.W.dot(&pi_a));
    }
    rates
  }

  /// Compute the transition probability matrices for all sites at time t.
  ///
  /// Returns P_a(t) = v_a * diag(exp(eigvals_a * mu_a * t)) * v_inv_a for each site a.
  /// Shape: [n_states, n_states, seq_len] where P[:,:,a] is the column-stochastic
  /// transition matrix for site a.
  ///
  /// Uses interpolation when available and t is within the interpolation range.
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

  /// Compute exp(Q*t) directly from eigendecomposition (no interpolation).
  ///
  /// P_a(t)[i,k] = sum_j v_a[i,j] * exp(eigvals_a[j] * mu_a * t) * v_inv_a[j,k]
  ///
  /// This is the computational bottleneck avoided by interpolation. Each site
  /// requires O(n^2) work for the two matrix-vector products.
  pub fn expQt_raw(&self, t: f64) -> Array3<f64> {
    let n = self.eigvals.nrows();
    let mut result = Array3::zeros((n, n, self.seq_len));

    for a in 0..self.seq_len {
      // Column scaling: v_a * diag(exp_lambda) = v_a[:, j] * exp_lambda[j]
      // avoids constructing an explicit diagonal matrix
      let e_lambda_t: Array1<f64> = (&self.eigvals.column(a) * self.mu[a] * t).mapv(f64::exp);
      let v_a = self.v.slice(s![.., .., a]);
      let v_inv_a = self.v_inv.slice(s![.., .., a]);
      let scaled_v = &v_a * &e_lambda_t;
      let p_a = scaled_v.dot(&v_inv_a);
      result.slice_mut(s![.., .., a]).assign(&clamp_min(&p_a, 0.0));
    }

    result
  }

  /// Propagate a sequence profile backward in time (child -> parent).
  ///
  /// For each site a: result[a, j] = sum_i profile[a, i] * P_a(t)[i, j]
  ///
  /// This is the Felsenstein pruning step: given likelihoods at a child node,
  /// compute the parent's partial likelihood contribution from this child.
  ///
  /// # Arguments
  ///
  /// * `profile` - Child state likelihoods. Shape: [seq_len, n_states].
  /// * `t` - Branch length (must be non-negative).
  /// * `return_log` - If true, return log-likelihoods.
  ///
  /// # Returns
  ///
  /// Parent partial likelihoods. Shape: [seq_len, n_states].
  pub fn propagate_profile(&self, profile: &Array2<f64>, t: f64, return_log: bool) -> Result<Array2<f64>, Report> {
    let qt = self.expQt(t)?;
    let mut result = Array2::zeros(profile.dim());

    // For each site a: result[a,:] = profile[a,:] @ Qt[:,:,a]
    for a in 0..self.seq_len {
      let qt_a = qt.slice(s![.., .., a]);
      result.row_mut(a).assign(&profile.row(a).dot(&qt_a));
    }

    if return_log {
      result.mapv_inplace(f64::ln);
    }
    Ok(result)
  }

  /// Evolve a sequence profile forward in time (parent -> child).
  ///
  /// For each site a: result[a, i] = sum_j profile[a, j] * P_a(t)[i, j]
  ///
  /// Given the parent's state distribution, compute the child's expected distribution
  /// after time t of evolution.
  ///
  /// # Arguments
  ///
  /// * `profile` - Parent state probabilities. Shape: [seq_len, n_states].
  /// * `t` - Branch length (must be non-negative).
  /// * `return_log` - If true, return log-probabilities.
  ///
  /// # Returns
  ///
  /// Child state probabilities. Shape: [seq_len, n_states].
  pub fn evolve(&self, profile: &Array2<f64>, t: f64, return_log: bool) -> Result<Array2<f64>, Report> {
    let qt = self.expQt(t)?;
    let mut result = Array2::zeros(profile.dim());

    // For each site a: result[a,:] = profile[a,:] @ Qt[:,:,a]^T
    for a in 0..self.seq_len {
      let qt_a = qt.slice(s![.., .., a]);
      result.row_mut(a).assign(&profile.row(a).dot(&qt_a.t()));
    }

    if return_log {
      result.mapv_inplace(f64::ln);
    }
    Ok(result)
  }

  /// Construct the per-site rate matrices Q_a in column-stochastic form for display.
  ///
  /// Returns Q_a where Q_a[i,j] = W[i,j] * pi_a[i] for i != j, and columns sum to 0.
  /// Shape: [n_states, n_states, seq_len].
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

  /// Build or rebuild the interpolation table.
  ///
  /// Pre-computes exp(Qt) on a non-uniform grid of t values, denser near t=0
  /// where the matrix changes most rapidly. The grid is scaled by the inverse
  /// of the mean substitution rate for numerical stability.
  fn build_interpolator(&mut self) {
    let avg_rates = self.average_rate();
    let rate_scale = (avg_rates.sum() / self.seq_len as f64).max(1e-10);
    let inv_rate = 1.0 / rate_scale;

    // Non-uniform grid: dense near 0, progressively sparser
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

    // Pre-compute expQt at each grid point (grid t values are always non-negative)
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

/// Pre-computed exp(Qt) interpolation table for efficient evaluation.
///
/// Stores exp(Qt) at discrete time grid points and uses linear interpolation
/// for arbitrary t values. This avoids repeated eigendecomposition during
/// tree traversal, reducing per-branch cost from O(L * n^3) to O(L * n^2).
#[derive(Clone, Debug)]
struct ExpQtInterpolator {
  /// Non-uniform time grid points (sorted ascending). Shape: [n_t].
  t_grid: Array1<f64>,
  /// Pre-computed matrices. Shape: [n_t, n_states, n_states, seq_len].
  data: Array4<f64>,
  /// Mean substitution rate, used to determine interpolation range.
  rate_scale: f64,
}

impl ExpQtInterpolator {
  /// Maximum scaled time value for which interpolation is used.
  /// Beyond this, fall back to direct computation.
  const MAX_INTERP_RANGE: f64 = 10.0;

  /// Linearly interpolate exp(Qt) at time t.
  ///
  /// Finds the bracketing grid interval and computes:
  ///   result = (1 - alpha) * data[left] + alpha * data[left + 1]
  fn interpolate(&self, t: f64) -> Array3<f64> {
    let n = self.t_grid.len();
    debug_assert!(n >= 2);

    // Clamp to grid range (extrapolation not meaningful for probability matrices)
    let t = t.clamp(self.t_grid[0], self.t_grid[n - 1]);

    // Linear scan for bracketing interval (61-element grid, O(n) is fine)
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

/// Generate linearly spaced values (inclusive of both endpoints).
fn linspace(start: f64, end: f64, n: usize) -> Vec<f64> {
  if n <= 1 {
    return vec![start];
  }
  let step = (end - start) / (n - 1) as f64;
  (0..n).map(|i| start + i as f64 * step).collect()
}
