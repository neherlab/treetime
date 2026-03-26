use crate::gtr::gtr_site_specific::GTRSiteSpecific;
use crate::gtr::infer_gtr::site_specific::MutationCountsSiteSpecific;
use ndarray::prelude::*;
use ndarray::Array3;

/// Simulate per-site mutation counts from a known site-specific GTR model.
///
/// For each site a, generates synthetic n_ija and T_ia consistent with the
/// model's W, pi_a, and mu_a: n_ija[i,j,a] = pi_a[i] * W[i,j] * T_ia[j,a] * mu_a
pub fn simulate_counts(gtr: &GTRSiteSpecific, total_time: f64) -> MutationCountsSiteSpecific {
  let n = gtr.pi.nrows();
  let seq_len = gtr.seq_len;

  let mut T_ia = Array2::zeros((n, seq_len));
  for a in 0..seq_len {
    for i in 0..n {
      T_ia[[i, a]] = total_time * gtr.pi[[i, a]];
    }
  }

  let mut n_ija = Array3::zeros((n, n, seq_len));
  for a in 0..seq_len {
    for i in 0..n {
      for j in 0..n {
        if i != j {
          n_ija[[i, j, a]] = gtr.pi[[i, a]] * gtr.W[[i, j]] * T_ia[[j, a]] * gtr.mu[a];
        }
      }
    }
  }

  let root_state = gtr.pi.clone();
  MutationCountsSiteSpecific {
    n_ija,
    T_ia,
    root_state,
  }
}

/// Convert nested JSON array [[r0c0, r0c1, ...], [r1c0, ...]] to Array2.
pub fn value_to_array2(value: &serde_json::Value) -> Array2<f64> {
  let rows: Vec<Vec<f64>> = serde_json::from_value(value.clone()).unwrap();
  let nrows = rows.len();
  let ncols = rows[0].len();
  let flat: Vec<f64> = rows.into_iter().flatten().collect();
  Array2::from_shape_vec((nrows, ncols), flat).unwrap()
}

/// Convert nested JSON array [[[d000, d001, ...], ...], ...] to Array3.
pub fn value_to_array3(value: &serde_json::Value) -> Array3<f64> {
  let dim0: Vec<Vec<Vec<f64>>> = serde_json::from_value(value.clone()).unwrap();
  let d0 = dim0.len();
  let d1 = dim0[0].len();
  let d2 = dim0[0][0].len();
  let flat: Vec<f64> = dim0.into_iter().flatten().flatten().collect();
  Array3::from_shape_vec((d0, d1, d2), flat).unwrap()
}
