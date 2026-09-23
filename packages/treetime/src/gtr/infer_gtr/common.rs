use crate::gtr::gtr::avg_transition;
use eyre::Report;
use log::warn;
use ndarray::{Array1, Array2, Array3, ArrayView1, Axis};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use treetime_utils::array::ndarray::{is_max_above, outer};
use treetime_utils::array::serde::{array1_as_vec, array1_from_vec, array2_as_vec, array2_from_vec};

const TINY_NUMBER: f64 = 1e-12;

pub(crate) fn infer_gtr_impl(counts: &MutationCounts, options: &InferGtrOptions) -> Result<InferGtrResult, Report> {
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

    W = (&(&nij.view() + &nij.t() + 2.0 * &pc_mat) / mu)
      / (&outer(&pi, Ti)? + &outer(Ti, &pi)? + TINY_NUMBER + 2.0 * &pc_mat);
    W.diag_mut().fill(0.0);
    W /= avg_transition(&W, &pi)?;

    if fixed_pi.is_none() {
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

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MutationCounts {
  #[serde(serialize_with = "array2_as_vec", deserialize_with = "array2_from_vec")]
  pub nij: Array2<f64>,

  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub Ti: Array1<f64>,

  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub root_state: Array1<f64>,
}

#[derive(Clone, Debug, SmartDefault)]
pub struct InferGtrOptions {
  pub fixed_pi: Option<Array1<f64>>,

  #[default = 1.0]
  pub pc: f64,

  #[default = 1e-5]
  pub dp: f64,

  #[default = 40]
  pub max_iter: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct InferGtrResult {
  #[serde(serialize_with = "array2_as_vec", deserialize_with = "array2_from_vec")]
  pub W: Array2<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub pi: Array1<f64>,
  pub mu: f64,
}

pub(crate) fn distance(pi_old: &Array1<f64>, pi: &Array1<f64>) -> f64 {
  (pi_old - pi).mapv(|x| x * x).sum().sqrt()
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn is_profile_informative(profile: ArrayView1<'_, f64>, n_states: usize) -> bool {
  let uniform_threshold = 1.0 / n_states as f64 + 1e-10;
  is_max_above(&profile, uniform_threshold)
}

pub(crate) fn get_branch_mutation_matrix(
  msg_to_child: &Array2<f64>,
  msg_to_parent: &Array2<f64>,
  exp_qt: &Array2<f64>,
) -> Array3<f64> {
  let (n_sites, n_states) = msg_to_parent.dim();
  let mut result = Array3::zeros((n_sites, n_states, n_states));

  for a in 0..n_sites {
    let mut site_sum = 0.0;

    for i in 0..n_states {
      for j in 0..n_states {
        let val = msg_to_parent[[a, i]] * exp_qt[[i, j]] * msg_to_child[[a, j]];
        result[[a, i, j]] = val;
        site_sum += val;
      }
    }

    if site_sum > 0.0 {
      for i in 0..n_states {
        for j in 0..n_states {
          result[[a, i, j]] /= site_sum;
        }
      }
    }
  }

  result
}

pub(crate) fn accumulate_mutation_counts(
  mut_stack: &Array3<f64>,
  branch_length: f64,
  nij: &mut Array2<f64>,
  Ti: &mut Array1<f64>,
) {
  let (n_sites, n_states, _) = mut_stack.dim();

  for i in 0..n_states {
    for j in 0..n_states {
      let mut sum = 0.0;
      for a in 0..n_sites {
        sum += mut_stack[[a, i, j]];
      }
      nij[[i, j]] += sum;
    }
  }

  for k in 0..n_states {
    let mut parent_sum = 0.0;
    let mut child_sum = 0.0;
    for a in 0..n_sites {
      for i in 0..n_states {
        parent_sum += mut_stack[[a, i, k]];
      }
      for j in 0..n_states {
        child_sum += mut_stack[[a, k, j]];
      }
    }
    Ti[k] += 0.5 * branch_length * (parent_sum + child_sum);
  }
}
