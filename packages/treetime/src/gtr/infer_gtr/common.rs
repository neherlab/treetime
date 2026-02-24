use crate::gtr::gtr::avg_transition;
use eyre::Report;
use log::warn;
use ndarray::{Array1, Array2, Axis};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use treetime_utils::array::ndarray::outer;
use treetime_utils::array::serde::{array1_as_vec, array1_from_vec, array2_as_vec, array2_from_vec};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MutationCounts {
  /// NxN matrix where each entry represents the observed number of transitions from state i to state j.
  #[serde(serialize_with = "array2_as_vec", deserialize_with = "array2_from_vec")]
  pub nij: Array2<f64>,

  /// An N-vector representing the total time spent in each character state.
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub Ti: Array1<f64>,

  /// An N-vector representing the state counts at the root of the phylogenetic tree.
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub root_state: Array1<f64>,
}

#[derive(Clone, Debug, SmartDefault)]
pub struct InferGtrOptions {
  /// Optional fixed equilibrium state frequencies. If `None`, then frequencies are estimated.
  pub fixed_pi: Option<Array1<f64>>,

  /// Pseudo-counts to regularize zero observations in transition data.
  #[default = 1.0]
  pub pc: f64,

  /// Convergence criterion for the iterative optimization.
  #[default = 1e-5]
  pub dp: f64,

  /// Maximum number of iterations allowed for the optimization process.
  #[default = 40]
  pub max_iter: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct InferGtrResult {
  /// Substitution attempt matrix calculated during the GTR inference.
  #[serde(serialize_with = "array2_as_vec", deserialize_with = "array2_from_vec")]
  pub W: Array2<f64>,
  /// Estimated equilibrium state frequencies.
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub pi: Array1<f64>,
  /// Scale factor for the rates matrix.
  pub mu: f64,
}

/// Small constant to prevent division by zero, matching Python's ttconf.TINY_NUMBER.
const TINY_NUMBER: f64 = 1e-12;

/// Infer a GTR model by specifying the number of transitions and time spent in each character. The basic equation
/// that is being solved is
///
/// $$ n_{ij} = pi_i W_{ij} T_j $$
///
/// where:
///
///  * $n_{ij}$ are the transitions
///  * $pi_i$ are the equilibrium state frequencies
///  * $W_{ij}$ is the "substitution attempt matrix"
///  * $T_i$ is the time on the tree spent in character state $i$.
///
/// To regularize the process, we add pseudo-counts and also need to account for the fact that the root of the tree is
/// in a particular state. the modified equation is
///
/// $$ n_{ij} + pc = pi_i W_{ij} (T_j+pc+root\_state) $$
pub fn infer_gtr_impl(counts: &MutationCounts, options: &InferGtrOptions) -> Result<InferGtrResult, Report> {
  let MutationCounts { nij, Ti, root_state } = counts;
  let InferGtrOptions {
    fixed_pi,
    pc,
    dp,
    max_iter,
  } = options;

  let N = Ti.len();

  let pc_mat = {
    let mut pc_mat = Array2::from_elem((N, N), *pc);
    pc_mat.diag_mut().fill(0.0);
    pc_mat
  };

  let nij = {
    let mut nij = nij.clone();
    nij.diag_mut().fill(0.0);
    nij
  };

  let mut pi_old = Array1::zeros(N);
  let mut pi = fixed_pi.clone().unwrap_or_else(|| Array1::ones(N));
  pi /= pi.sum();
  let mut W = Array2::ones((N, N));
  let mut mu = (nij.sum() + pc) / (Ti.sum() + pc);

  for _ in 0..*max_iter {
    let dist = distance(&pi_old, &pi);

    if dist < *dp {
      break;
    }

    pi_old.assign(&pi);

    // W calculation matching Python: includes TINY_NUMBER for numerical stability
    W = (&(&nij.view() + &nij.t() + 2.0 * &pc_mat) / mu)
      / (&outer(&pi, Ti)? + &outer(Ti, &pi)? + TINY_NUMBER + 2.0 * &pc_mat);
    W.diag_mut().fill(0.0);
    W /= avg_transition(&W, &pi)?;

    if fixed_pi.is_none() {
      // pi calculation matching Python: includes TINY_NUMBER for numerical stability
      pi = (&nij.sum_axis(Axis(1)) + &pc_mat.sum_axis(Axis(1)) + root_state)
        / (TINY_NUMBER + mu * W.dot(Ti) + root_state.sum() + pc_mat.sum_axis(Axis(1)));
      pi /= pi.sum();
      mu = (nij.sum() + pc) / (pi.dot(&W.dot(Ti)) + pc);
    } else {
      mu = (nij.sum() + pc) / (pi.dot(&W.dot(&pi)) * Ti.sum() + pc);
    }
  }

  if distance(&pi_old, &pi) > *dp {
    warn!("When inferring GTR parameters: The iterative scheme has not converged.");
  } else if (pi.sum() - 1.0).abs() > *dp {
    warn!("When inferring GTR parameters: Proper normalization was not reached.");
  }
  Ok(InferGtrResult { W, pi, mu })
}

pub fn distance(pi_old: &Array1<f64>, pi: &Array1<f64>) -> f64 {
  (pi_old - pi).mapv(|x| x * x).sum().sqrt()
}
