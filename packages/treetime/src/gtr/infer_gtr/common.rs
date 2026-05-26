use crate::gtr::gtr::avg_transition;
use eyre::Report;
use log::warn;
use ndarray::{Array1, Array2, Array3, ArrayView1, Axis};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use treetime_utils::array::ndarray::{is_max_above, outer};
use treetime_utils::array::serde::{array1_as_vec, array1_from_vec, array2_as_vec, array2_from_vec};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MutationCounts {
  /// Expected substitution counts accumulated across all edges and sites.
  ///
  /// `nij[i, j]`: expected count of substitutions from parent state `j` to child state `i`.
  ///
  /// - Dense: fractional posterior mass from the branch joint distribution
  ///   `P(child=i, parent=j | site, data)`, summed over sites and edges.
  /// - Sparse: integer count of observed mutations (from Fitch reconstruction),
  ///   one per substitution event per edge.
  ///
  /// Diagonal is zeroed after accumulation (no-change events excluded).
  #[serde(serialize_with = "array2_as_vec", deserialize_with = "array2_from_vec")]
  pub nij: Array2<f64>,

  /// Total evolutionary time spent in each character state across the tree.
  ///
  /// `Ti[k]`: total time in state `k`, summed over all edges and sites.
  ///
  /// - Dense: midpoint approximation per edge:
  ///   `Ti[k] += 0.5 * branch_length * (P(parent=k) + P(child=k))`
  ///   where marginals come from summing the branch joint distribution.
  /// - Sparse: `branch_length * child_composition_count(k)`,
  ///   adjusted by `-0.5*BL` at mutation sites (child state loses half)
  ///   and `+0.5*BL` (parent state gains half).
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub Ti: Array1<f64>,

  /// State composition at the tree root, used as a prior on equilibrium frequencies.
  ///
  /// `root_state[k]`: count of positions in state `k` at the root.
  ///
  /// - Dense: argmax of marginal profile per site, counting most-likely states.
  /// - Sparse: composition counts from Fitch consensus sequence.
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
/// in a particular state. The modified equation is
///
/// $$ n_{ij} + pc = pi_i W_{ij} (T_j+pc+root\_state) $$
///
/// with two regularization terms:
///
///  * $pc$ is a scalar pseudo-count that regularizes zero observations in transition data
///  * $root\_state$ is a vector of character-state counts at the root, acting as a prior on equilibrium frequencies
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

/// Whether a profile is peaked above the uniform baseline (carries phylogenetic signal).
pub fn is_profile_informative(profile: ArrayView1<'_, f64>, n_states: usize) -> bool {
  let uniform_threshold = 1.0 / n_states as f64 + 1e-10;
  is_max_above(&profile, uniform_threshold)
}

/// Compute branch joint distribution for GTR inference.
///
/// Given edge messages and transition matrix, computes the posterior probability
/// P(child=i, parent=j | site a) for each site.
///
/// Arguments:
/// - `msg_to_child`: Outgroup likelihood at child, indexed by parent state (L x n)
/// - `msg_to_parent`: Subtree likelihood at child, indexed by child state (L x n)
/// - `exp_qt`: Transition matrix P(child=i | parent=j) (n x n)
///
/// Returns: Array3 of shape (L, n, n) where result[a, i, j] = P(child=i, parent=j | site a)
pub fn get_branch_mutation_matrix(
  msg_to_child: &Array2<f64>,
  msg_to_parent: &Array2<f64>,
  exp_qt: &Array2<f64>,
) -> Array3<f64> {
  let (n_sites, n_states) = msg_to_parent.dim();
  let mut result = Array3::zeros((n_sites, n_states, n_states));

  for a in 0..n_sites {
    let mut site_sum = 0.0;

    // Compute unnormalized joint: pc[a,i] * exp_qt[i,j] * pp[a,j]
    for i in 0..n_states {
      for j in 0..n_states {
        let val = msg_to_parent[[a, i]] * exp_qt[[i, j]] * msg_to_child[[a, j]];
        result[[a, i, j]] = val;
        site_sum += val;
      }
    }

    // Normalize to get posterior probabilities
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

/// Accumulate mutation counts and time-in-state from branch joint distribution.
///
/// Updates:
/// - `nij`: adds expected substitution counts summed over sites
/// - `Ti`: adds time spent in each state using midpoint approximation
///
/// The midpoint formula: Ti += 0.5 * branch_length * (P(parent=i) + P(child=i))
/// where marginals are computed by summing mut_stack over child/parent states.
pub fn accumulate_mutation_counts(
  mut_stack: &Array3<f64>,
  branch_length: f64,
  nij: &mut Array2<f64>,
  Ti: &mut Array1<f64>,
) {
  let (n_sites, n_states, _) = mut_stack.dim();

  // nij += sum over sites: mut_stack.sum(axis=0)
  for i in 0..n_states {
    for j in 0..n_states {
      let mut sum = 0.0;
      for a in 0..n_sites {
        sum += mut_stack[[a, i, j]];
      }
      nij[[i, j]] += sum;
    }
  }

  // Compute marginals and accumulate Ti
  // parent_marginal[a, j] = sum_i mut_stack[a, i, j] (sum over child states)
  // child_marginal[a, i] = sum_j mut_stack[a, i, j] (sum over parent states)
  // Ti[k] += 0.5 * branch_length * (sum_a parent_marginal[a,k] + sum_a child_marginal[a,k])
  for k in 0..n_states {
    let mut parent_sum = 0.0;
    let mut child_sum = 0.0;
    for a in 0..n_sites {
      // parent_marginal: sum over child states i
      for i in 0..n_states {
        parent_sum += mut_stack[[a, i, k]];
      }
      // child_marginal: sum over parent states j
      for j in 0..n_states {
        child_sum += mut_stack[[a, k, j]];
      }
    }
    Ti[k] += 0.5 * branch_length * (parent_sum + child_sum);
  }
}
