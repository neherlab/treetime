use crate::gtr::gtr_site_specific::{GTRSiteSpecific, GTRSiteSpecificParams};
use eyre::Report;
use log::warn;
use ndarray::Array3;
use ndarray::prelude::*;
use smart_default::SmartDefault;

/// Per-site mutation counts for site-specific GTR inference.
///
/// Unlike `MutationCounts` which sums over sites, this preserves per-site
/// statistics needed to estimate site-specific equilibrium frequencies.
#[derive(Clone, Debug)]
pub struct MutationCountsSiteSpecific {
  /// Per-site substitution counts. Shape: [n_states, n_states, seq_len].
  /// `n_ija[i, j, a]`: expected count of substitutions from parent state j
  /// to child state i at site a.
  pub n_ija: Array3<f64>,

  /// Per-site time-in-state. Shape: [n_states, seq_len].
  /// `T_ia[i, a]`: total evolutionary time spent in state i at site a.
  pub T_ia: Array2<f64>,

  /// Per-site root state distribution. Shape: [n_states, seq_len].
  /// `root_state[i, a]`: probability that site a is in state i at the root.
  pub root_state: Array2<f64>,
}

/// Options for site-specific GTR inference.
#[derive(Clone, Debug, SmartDefault)]
pub struct InferGtrSiteSpecificOptions {
  /// Number of character states.
  pub n_states: usize,

  /// Pseudo-counts for regularization.
  #[default = 1.0]
  pub pc: f64,

  /// Lower bound on gap state frequency (prevents artifacts from near-zero gap rates).
  #[default = 0.01]
  pub gap_limit: f64,

  /// Index of the gap character in the alphabet (None if no gap character).
  pub gap_index: Option<usize>,

  /// Maximum number of iterations.
  #[default = 30]
  pub max_iter: usize,

  /// Convergence threshold (L2 norm of pi change).
  #[default = 1e-5]
  pub dp: f64,
}

/// Result of site-specific GTR inference.
#[derive(Clone, Debug)]
pub struct InferGtrSiteSpecificResult {
  /// Shared exchangeability matrix. Shape: [n_states, n_states].
  pub W: Array2<f64>,
  /// Per-site equilibrium frequencies. Shape: [n_states, seq_len].
  pub pi: Array2<f64>,
  /// Per-site substitution rates. Shape: [seq_len].
  pub mu: Array1<f64>,
}

/// Infer site-specific GTR model parameters from per-site mutation statistics.
///
/// Solves the equation n_ij,a = pi_i,a * W_ij * T_j,a * mu_a iteratively:
///
/// 1. Estimate W from aggregated (cross-site) mutation counts and expected counts.
/// 2. Estimate per-site pi from per-site mutation marginals and time-in-state.
/// 3. Estimate per-site mu from total per-site mutation rate.
/// 4. Repeat until convergence.
///
/// The shared W captures relative exchangeability between states (e.g. transition/transversion
/// bias). Per-site pi captures compositional heterogeneity across the alignment (e.g. conserved
/// positions vs variable positions). Per-site mu captures rate heterogeneity.
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
    ..
  } = options;

  let q = *n_states;
  let seq_len = n_ija.shape()[2];

  // Zero the diagonal of n_ija (self-transitions are not mutations)
  let n_ija = {
    let mut n_ija = n_ija.clone();
    for a in 0..seq_len {
      for i in 0..q {
        n_ija[[i, i, a]] = 0.0;
      }
    }
    n_ija
  };

  // n_ij: site-summed mutation counts [n_states, n_states]
  let n_ij = n_ija.sum_axis(Axis(2));

  // m_ia: marginal mutation counts per site + root + pseudocounts [n_states, seq_len]
  // sum over parent states (axis=1) gives total mutations into each child state per site
  let m_ia = n_ija.sum_axis(Axis(1)) + root_state + *pc;

  // n_a: total mutations per site + pseudocounts [seq_len]
  let n_a = n_ija.sum_axis(Axis(1)).sum_axis(Axis(0)) + *pc;

  // Lambda: normalization constant per site (total root probability mass + q * pc)
  let Lambda = root_state.sum_axis(Axis(0)) + q as f64 * *pc;

  // Initialize iterates
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

    // S_ij = sum_a mu_a * p_ia[:,a] outer T_ia[:,a]
    // S_ij[i,j] = sum_a mu_a * p_ia[i,a] * T_ia[j,a]
    let S_ij = einsum_mu_pi_T(&mu_a, &p_ia, T_ia);

    // W_ij = (n_ij + n_ij^T + pc) / (S_ij + S_ij^T + pc)
    W_ij = (&n_ij + &n_ij.t() + *pc) / (&S_ij + &S_ij.t() + *pc);
    W_ij.diag_mut().fill(0.0);

    // Normalize W by average rate (using mean pi across sites)
    let avg_pi: Array1<f64> = p_ia.sum_axis(Axis(1)) / seq_len as f64;
    let average_rate = avg_pi.dot(&W_ij.dot(&avg_pi));
    if average_rate > 0.0 {
      W_ij /= average_rate;
      mu_a *= average_rate;
    }

    // Update per-site frequencies: p_ia = m_ia / (mu_a * W @ T_ia + Lambda)
    let W_T_ia = W_ij.dot(T_ia);
    p_ia = &m_ia / (&W_T_ia * &mu_a + &Lambda);
    // Normalize columns to sum to 1
    let col_sums = p_ia.sum_axis(Axis(0));
    p_ia /= &col_sums;

    // Update per-site rates: mu_a = n_a / (pc + einsum('ia,ij,ja->a', p_ia, W_ij, T_ia))
    let denominator = einsum_pi_W_T(&p_ia, &W_ij, T_ia) + *pc;
    mu_a = &n_a / &denominator;
  }

  if l2_norm_diff(&p_ia_old, &p_ia) > *dp {
    warn!("Site-specific GTR inference: maximum iterations reached without convergence.");
  }

  // Enforce gap frequency lower bound
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

/// Construct a `GTRSiteSpecific` model from inference results.
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

/// S_ij = sum_a mu_a * p_ia[i,a] * T_ia[j,a]
///
/// Equivalent to numpy: `np.einsum('a,ia,ja', mu_a, p_ia, T_ia)`
fn einsum_mu_pi_T(mu_a: &Array1<f64>, p_ia: &Array2<f64>, T_ia: &Array2<f64>) -> Array2<f64> {
  let (q, seq_len) = p_ia.dim();
  let mut result = Array2::zeros((q, q));
  for a in 0..seq_len {
    let mu = mu_a[a];
    let p_col = p_ia.column(a);
    let t_col = T_ia.column(a);
    // result += mu * outer(p_col, t_col)
    for i in 0..q {
      for j in 0..q {
        result[[i, j]] += mu * p_col[i] * t_col[j];
      }
    }
  }
  result
}

/// result_a = sum_{i,j} p_ia[i,a] * W_ij[i,j] * T_ia[j,a]
///
/// Equivalent to numpy: `np.einsum('ia,ij,ja->a', p_ia, W_ij, T_ia)`
fn einsum_pi_W_T(p_ia: &Array2<f64>, W_ij: &Array2<f64>, T_ia: &Array2<f64>) -> Array1<f64> {
  let (q, seq_len) = p_ia.dim();
  let mut result = Array1::zeros(seq_len);
  // Precompute W @ T_ia to avoid redundant inner loop
  let W_T = W_ij.dot(T_ia); // [q, seq_len]
  for a in 0..seq_len {
    result[a] = p_ia.column(a).dot(&W_T.column(a));
  }
  result
}

/// L2 norm of the difference between two 2D arrays.
fn l2_norm_diff(a: &Array2<f64>, b: &Array2<f64>) -> f64 {
  (a - b).mapv(|x| x * x).sum().sqrt()
}
