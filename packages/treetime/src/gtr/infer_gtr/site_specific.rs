use crate::gtr::gtr_site_specific::{GTRSiteSpecific, GTRSiteSpecificParams};
use eyre::Report;
use log::warn;
use ndarray::Array3;
use ndarray::prelude::*;
use smart_default::SmartDefault;

#[derive(Clone, Debug)]
pub struct MutationCountsSiteSpecific {
  pub n_ija: Array3<f64>,

  pub T_ia: Array2<f64>,

  pub root_state: Array2<f64>,
}

#[derive(Clone, Debug, SmartDefault)]
pub struct InferGtrSiteSpecificOptions {
  pub n_states: usize,

  #[default = 1.0]
  pub pc: f64,

  #[default = 0.01]
  pub gap_limit: f64,

  pub gap_index: Option<usize>,

  #[default = 30]
  pub max_iter: usize,

  #[default = 1e-5]
  pub dp: f64,
}

#[derive(Clone, Debug)]
pub struct InferGtrSiteSpecificResult {
  pub W: Array2<f64>,
  pub pi: Array2<f64>,
  pub mu: Array1<f64>,
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub fn infer_gtr_site_specific_impl(
  counts: &MutationCountsSiteSpecific,
  options: &InferGtrSiteSpecificOptions,
) -> Result<InferGtrSiteSpecificResult, Report> {
  let MutationCountsSiteSpecific {
    n_ija,
    T_ia,
    root_state,
  } = counts;
  let InferGtrSiteSpecificOptions {
    n_states,
    pc,
    gap_limit,
    gap_index,
    max_iter,
    dp,
  } = options;

  let q = *n_states;
  let seq_len = n_ija.shape()[2];

  let n_ija = {
    let mut n_ija = n_ija.clone();
    for a in 0..seq_len {
      for i in 0..q {
        n_ija[[i, i, a]] = 0.0;
      }
    }
    n_ija
  };

  let n_ij = n_ija.sum_axis(Axis(2));

  let m_ia = n_ija.sum_axis(Axis(1)) + root_state + *pc;

  let n_a = n_ija.sum_axis(Axis(1)).sum_axis(Axis(0)) + *pc;

  let Lambda = root_state.sum_axis(Axis(0)) + q as f64 * *pc;

  let mut p_ia_old = Array2::zeros((q, seq_len));
  let mut p_ia = Array2::from_elem((q, seq_len), 1.0 / q as f64);
  let mut mu_a = Array1::ones(seq_len);
  let mut W_ij = {
    let mut w = Array2::ones((q, q));
    w.diag_mut().fill(0.0);
    w
  };

  for iter in 0..*max_iter {
    let dist = l2_norm_diff(&p_ia_old, &p_ia);
    if iter > 0 && dist < *dp {
      break;
    }
    p_ia_old.assign(&p_ia);

    let S_ij = einsum_mu_pi_T(&mu_a, &p_ia, T_ia);

    W_ij = (&n_ij + &n_ij.t() + *pc) / (&S_ij + &S_ij.t() + *pc);
    W_ij.diag_mut().fill(0.0);

    let avg_pi: Array1<f64> = p_ia.sum_axis(Axis(1)) / seq_len as f64;
    let average_rate = avg_pi.dot(&W_ij.dot(&avg_pi));
    if average_rate > 0.0 {
      W_ij /= average_rate;
      mu_a *= average_rate;
    }

    let W_T_ia = W_ij.dot(T_ia);
    p_ia = &m_ia / (&W_T_ia * &mu_a + &Lambda);
    let col_sums = p_ia.sum_axis(Axis(0));
    p_ia /= &col_sums;

    let denominator = einsum_pi_W_T(&p_ia, &W_ij, T_ia) + *pc;
    mu_a = &n_a / &denominator;
  }

  if l2_norm_diff(&p_ia_old, &p_ia) > *dp {
    warn!("Site-specific GTR inference: maximum iterations reached without convergence.");
  }

  if let Some(gap_idx) = gap_index {
    for a in 0..seq_len {
      if p_ia[[*gap_idx, a]] < *gap_limit {
        p_ia[[*gap_idx, a]] = *gap_limit;
        let col_sum = p_ia.column(a).sum();
        p_ia.column_mut(a).mapv_inplace(|v| v / col_sum);
      }
    }
  }

  Ok(InferGtrSiteSpecificResult {
    W: W_ij,
    pi: p_ia,
    mu: mu_a,
  })
}

pub fn build_gtr_site_specific(
  result: &InferGtrSiteSpecificResult,
  n_states: usize,
  approximate: bool,
) -> Result<GTRSiteSpecific, Report> {
  let seq_len = result.pi.ncols();
  GTRSiteSpecific::new(GTRSiteSpecificParams {
    n_states,
    seq_len,
    mu: result.mu.clone(),
    W: Some(result.W.clone()),
    pi: result.pi.clone(),
    approximate,
  })
}

pub fn einsum_mu_pi_T(mu_a: &Array1<f64>, p_ia: &Array2<f64>, T_ia: &Array2<f64>) -> Array2<f64> {
  let (q, seq_len) = p_ia.dim();
  let mut result = Array2::zeros((q, q));
  for a in 0..seq_len {
    let mu = mu_a[a];
    let p_col = p_ia.column(a);
    let t_col = T_ia.column(a);
    for i in 0..q {
      for j in 0..q {
        result[[i, j]] += mu * p_col[i] * t_col[j];
      }
    }
  }
  result
}

pub fn einsum_pi_W_T(p_ia: &Array2<f64>, W_ij: &Array2<f64>, T_ia: &Array2<f64>) -> Array1<f64> {
  let (q, seq_len) = p_ia.dim();
  let mut result = Array1::zeros(seq_len);
  let W_T = W_ij.dot(T_ia);
  for a in 0..seq_len {
    result[a] = p_ia.column(a).dot(&W_T.column(a));
  }
  result
}

pub fn l2_norm_diff(a: &Array2<f64>, b: &Array2<f64>) -> f64 {
  (a - b).mapv(|x| x * x).sum().sqrt()
}
