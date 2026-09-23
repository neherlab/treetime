use crate::make_error;
use bon::bon;
use eyre::Report;
use ndarray::prelude::*;
use ndarray_linalg::Eigh;
use ndarray_linalg::UPLO::Lower;
use num_traits::abs;
use serde::Serialize;
use treetime_utils::array::ndarray::{clamp_min, outer};
use treetime_utils::array::serde::{array1_as_vec, array2_as_vec, option_array1_as_vec};

#[derive(Clone, Debug, Serialize)]
pub struct GTR {
  pub debug: bool,
  pub average_rate: f64,
  pub mu: f64,
  #[serde(serialize_with = "array2_as_vec")]
  pub W: Array2<f64>,
  #[serde(serialize_with = "array1_as_vec")]
  pub pi: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec")]
  pub eigvals: Array1<f64>,
  #[serde(serialize_with = "array2_as_vec")]
  pub v: Array2<f64>,
  #[serde(serialize_with = "array2_as_vec")]
  pub v_inv: Array2<f64>,
  #[serde(serialize_with = "option_array1_as_vec")]
  pub site_rates: Option<Array1<f64>>,
  pub unimodal_branch_likelihood: bool,
}

#[bon]
impl GTR {
  #[builder]
  pub fn new(n_states: usize, mu: f64, W: Option<Array2<f64>>, pi: Array1<f64>) -> Result<Self, Report> {
    let n = n_states;

    if pi.shape() != [n] {
      return make_error!(
        "Length of equilibrium frequency vector (pi) is {}, expected {n}",
        pi.len()
      );
    }

    if let Some(W) = &W {
      if W.shape() != [n, n] {
        return make_error!(
          "Dimensions of substitution matrix (W) are {:?}, expected [{n}, {n}]",
          W.shape()
        );
      }
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

    let pi = {
      let pi_sum = pi.sum();
      pi / pi_sum
    };

    let average_rate = avg_transition(&W, &pi)?;
    let mu = mu * average_rate;
    let W = W / average_rate;

    let (eigvals, v, v_inv) = eig_single_site(&W, pi.view())?;

    let unimodal_branch_likelihood = n == 2;

    Ok(Self {
      debug: false,
      average_rate,
      mu,
      W,
      pi,
      eigvals,
      v,
      v_inv,
      site_rates: None,
      unimodal_branch_likelihood,
    })
  }

  pub const fn average_rate(&self) -> f64 {
    self.average_rate
  }

  pub(crate) fn has_site_rates(&self) -> bool {
    self.site_rates.is_some()
  }

  pub fn set_site_rates(&mut self, rates: Array1<f64>) {
    self.site_rates = Some(rates);
  }

  pub fn clear_site_rates(&mut self) {
    self.site_rates = None;
  }

  pub(crate) fn expQt_with_rate(&self, t: f64, rate: f64) -> Array2<f64> {
    let exp_lt = self.exp_lt_scaled(t, rate);
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
        let scaled_eigvals = &self.eigvals * (self.mu * t);
        let exp_lt = (&rates.view().insert_axis(Axis(1)) * &scaled_eigvals.view().insert_axis(Axis(0))).mapv(f64::exp);
        let scaled = &transformed * &exp_lt;
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
        let scaled_eigvals = &self.eigvals * (self.mu * t);
        let exp_lt = (&rates.view().insert_axis(Axis(1)) * &scaled_eigvals.view().insert_axis(Axis(0))).mapv(f64::exp);
        let scaled = &transformed * &exp_lt;
        clamp_min(&scaled.dot(&self.v_inv), 0.0)
      },
    };
    if return_log { res.mapv(f64::ln) } else { res }
  }

  pub(crate) fn expQt(&self, t: f64) -> Array2<f64> {
    let exp_lt = self.exp_lt(t);
    let scaled_v_inv = &self.v_inv * &exp_lt.view().insert_axis(Axis(1));
    let Qt = self.v.dot(&scaled_v_inv);

    clamp_min(&Qt, 0.0)
  }

  fn exp_lt(&self, t: f64) -> Array1<f64> {
    (self.mu * t * &self.eigvals).mapv(f64::exp)
  }

  fn exp_lt_scaled(&self, t: f64, rate: f64) -> Array1<f64> {
    (self.mu * rate * t * &self.eigvals).mapv(f64::exp)
  }

  pub fn exp_eigvals_branch_length(&self, branch_length: f64) -> Array1<f64> {
    (&self.eigvals * branch_length).mapv(f64::exp)
  }

  pub fn Q(&self) -> Array2<f64> {
    let mut Q = (&self.W * &self.pi).t().to_owned();
    let diag = -Q.sum_axis(Axis(0));
    Q.diag_mut().assign(&diag);
    Q
  }
}

pub(crate) fn avg_transition(W: &Array2<f64>, pi: &Array1<f64>) -> Result<f64, Report> {
  Ok(pi.dot(W).dot(pi))
}

pub(super) fn eig_single_site(
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
