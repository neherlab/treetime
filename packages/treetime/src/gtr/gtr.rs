use crate::make_error;
use bon::bon;
use eyre::Report;
use getset::{CopyGetters, Getters};
use ndarray::prelude::*;
use ndarray_linalg::Eigh;
use ndarray_linalg::UPLO::Lower;
use num_traits::abs;
use treetime_utils::array::ndarray::{clamp_min, outer};

#[derive(Clone, Debug, Getters, CopyGetters)]
pub struct GTR {
  #[getset(get_copy = "pub")]
  mu: f64,
  #[getset(get = "pub")]
  W: Array2<f64>,
  #[getset(get = "pub")]
  pi: Array1<f64>,
  #[getset(get = "pub")]
  eigvals: Array1<f64>,
  #[getset(get = "pub")]
  v: Array2<f64>,
  #[getset(get = "pub")]
  v_inv: Array2<f64>,
  site_rates: Option<Array1<f64>>,
  unit_decay_rate: Option<f64>,
}

#[bon]
impl GTR {
  #[builder]
  pub fn new(
    n_states: usize,
    mu: f64,
    W: Option<Array2<f64>>,
    pi: Array1<f64>,
    rate_scale: RateScale,
    site_rates: Option<Array1<f64>>,
  ) -> Result<Self, Report> {
    let n = n_states;

    if pi.shape() != [n] {
      return make_error!(
        "Length of equilibrium frequency vector (pi) is {}, expected {n}",
        pi.len()
      );
    }

    if !pi.iter().all(|&p| p.is_finite() && p >= 0.0) || pi.sum() <= 0.0 {
      return make_error!("Equilibrium frequencies (pi) must be finite and non-negative with a positive sum, got {pi}");
    }

    if let Some(W) = &W {
      if W.shape() != [n, n] {
        return make_error!(
          "Dimensions of substitution matrix (W) are {:?}, expected [{n}, {n}]",
          W.shape()
        );
      }
      if !W.iter().all(|&w| w.is_finite() && w >= 0.0) {
        return make_error!("Substitution matrix (W) must contain finite non-negative entries, got {W}");
      }
    }

    validate_mu(mu)?;

    if let Some(site_rates) = &site_rates {
      validate_site_rates(site_rates)?;
    }

    let W = {
      let W = W.unwrap_or_else(|| {
        let mut W = Array2::<f64>::ones([n, n]);
        W.diag_mut().fill(0.0);
        W
      });
      let mut W = 0.5 * (&W.view() + &W.t());
      W.diag_mut().fill(0.0);
      W
    };

    let has_single_decay_rate = has_equal_positive_exchangeabilities(&W);

    let pi = {
      let pi_sum = pi.sum();
      pi / pi_sum
    };

    let average_rate = avg_transition(&W, &pi);
    if !(average_rate.is_finite() && average_rate > 0.0) {
      return make_error!(
        "Average substitution rate of the substitution matrix (W) under the equilibrium frequencies (pi) must be positive, got {average_rate}"
      );
    }
    let mu = match rate_scale {
      RateScale::FromMatrix => mu * average_rate,
      RateScale::Normalized => mu,
    };
    let W = W / average_rate;

    let unit_decay_rate = has_single_decay_rate.then(|| W[[0, 1]]);

    let (eigvals, v, v_inv) = eig_single_site(&W, pi.view())?;

    Ok(Self {
      mu,
      W,
      pi,
      eigvals,
      v,
      v_inv,
      site_rates,
      unit_decay_rate,
    })
  }

  pub fn site_rates(&self) -> Option<&Array1<f64>> {
    self.site_rates.as_ref()
  }

  pub fn certified_decay_rate(&self) -> Option<f64> {
    self.unit_decay_rate.map(|unit_decay_rate| unit_decay_rate * self.mu)
  }

  pub(crate) fn set_mu(&mut self, mu: f64) -> Result<(), Report> {
    validate_mu(mu)?;
    self.mu = mu;
    Ok(())
  }

  pub(crate) fn rate_eigvals(&self) -> Array1<f64> {
    self.rate_eigvals_at(1.0, 1.0)
  }

  pub(crate) fn expQt_with_rate(&self, t: f64, rate: f64) -> Array2<f64> {
    let exp_lt = self.exp_rate_eigvals(rate, t);
    let scaled_v_inv = &self.v_inv * &exp_lt.view().insert_axis(Axis(1));
    let Qt = self.v.dot(&scaled_v_inv);
    clamp_min(&Qt, 0.0)
  }

  pub(crate) fn evolve(&self, profile: &Array2<f64>, t: f64, return_log: bool) -> Array2<f64> {
    let res = match &self.site_rates {
      None => {
        let Qt = self.expQt(t);
        profile.dot(&Qt.t())
      },
      Some(rates) => {
        let transformed = profile.dot(&self.v_inv.t());
        let scaled = &transformed * &self.exp_site_rate_eigvals(rates, t);
        clamp_min(&scaled.dot(&self.v.t()), 0.0)
      },
    };
    if return_log { res.mapv(f64::ln) } else { res }
  }

  pub(crate) fn propagate_profile(&self, profile: &Array2<f64>, t: f64, return_log: bool) -> Array2<f64> {
    let res = match &self.site_rates {
      None => {
        let Qt = self.expQt(t);
        profile.dot(&Qt)
      },
      Some(rates) => {
        let transformed = profile.dot(&self.v);
        let scaled = &transformed * &self.exp_site_rate_eigvals(rates, t);
        clamp_min(&scaled.dot(&self.v_inv), 0.0)
      },
    };
    if return_log { res.mapv(f64::ln) } else { res }
  }

  pub(crate) fn expQt(&self, t: f64) -> Array2<f64> {
    self.expQt_with_rate(t, 1.0)
  }

  fn rate_eigvals_at(&self, rate: f64, t: f64) -> Array1<f64> {
    self.mu * rate * t * &self.eigvals
  }

  fn exp_rate_eigvals(&self, rate: f64, t: f64) -> Array1<f64> {
    self.rate_eigvals_at(rate, t).mapv_into(f64::exp)
  }

  fn exp_site_rate_eigvals(&self, rates: &Array1<f64>, t: f64) -> Array2<f64> {
    let rate_eigvals = self.rate_eigvals_at(1.0, t);
    (&rates.view().insert_axis(Axis(1)) * &rate_eigvals.view().insert_axis(Axis(0))).mapv_into(f64::exp)
  }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum RateScale {
  FromMatrix,
  Normalized,
}

pub(crate) fn avg_transition(W: &Array2<f64>, pi: &Array1<f64>) -> f64 {
  pi.dot(W).dot(pi)
}

fn validate_mu(mu: f64) -> Result<(), Report> {
  if mu.is_finite() && mu >= 0.0 {
    Ok(())
  } else {
    make_error!("Substitution rate (mu) must be finite and non-negative, got {mu}")
  }
}

fn validate_site_rates(site_rates: &Array1<f64>) -> Result<(), Report> {
  match site_rates.iter().position(|&rate| !(rate.is_finite() && rate >= 0.0)) {
    None => Ok(()),
    Some(site) => make_error!(
      "Site rates must be finite and non-negative, got {} at site {site}",
      site_rates[site]
    ),
  }
}

#[expect(
  clippy::float_cmp,
  reason = "the single decay rate certificate requires exactly equal exchangeabilities, without tolerance"
)]
fn has_equal_positive_exchangeabilities(W: &Array2<f64>) -> bool {
  let mut off_diagonal = W.indexed_iter().filter(|&((i, j), _)| i != j).map(|(_, &w)| w);
  off_diagonal
    .next()
    .is_some_and(|first| first > 0.0 && off_diagonal.all(|w| w == first))
}

fn eig_single_site(
  W: &Array2<f64>,
  pi: ArrayView1<'_, f64>,
) -> Result<(Array1<f64>, Array2<f64>, Array2<f64>), Report> {
  assert!(abs(W.diag().sum()) < 1e-10);

  let sqrt_pi: Array1<f64> = pi.mapv(f64::sqrt);
  let mut sym_Q: Array2<f64> = W * outer(&sqrt_pi, &sqrt_pi)?;

  let diag = -W.dot(&pi);
  sym_Q.diag_mut().assign(&diag);

  let (eigvals, eigvecs) = sym_Q.eigh(Lower)?;

  let tmp_v: Array2<f64> = eigvecs.t().to_owned() * sqrt_pi.to_owned();
  let one_norm: Array1<f64> = tmp_v.mapv(f64::abs).sum_axis(Axis(1));

  let v = tmp_v.t().to_owned() / &one_norm;
  let v_inv = (eigvecs * one_norm).t().to_owned() / sqrt_pi;

  Ok((eigvals, v, v_inv))
}
