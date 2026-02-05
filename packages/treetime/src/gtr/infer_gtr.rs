use crate::alphabet::alphabet::Alphabet;
use crate::graph::edge::{GraphEdgeKey, HasBranchLength};
use crate::graph::node::GraphNodeKey;
use crate::gtr::gtr::{GTR, GTRParams, avg_transition};
use crate::representation::graph_ancestral::GraphAncestral;
use crate::seq::composition::Composition;
use crate::seq::mutation::Sub;
use eyre::Report;
use log::warn;
use ndarray::{Array1, Array2, Axis};
use parking_lot::RwLock;
use smart_default::SmartDefault;
use std::sync::Arc;
use treetime_utils::ndarray::outer;

pub trait PartitionWithGtrInference {
  fn alphabet(&self) -> &Alphabet;
  fn get_seq_composition(&self, node_key: GraphNodeKey) -> &Composition;
  fn get_edge_substitutions(&self, edge_key: GraphEdgeKey, graph: &GraphAncestral) -> Vec<Sub>;
}

pub fn infer_gtr<P: PartitionWithGtrInference>(
  partition: &Arc<RwLock<P>>,
  graph: &GraphAncestral,
) -> Result<Result<GTR, Report>, Report> {
  let counts = get_mutation_counts(graph, partition)?;
  let InferGtrResult { W, pi, mu } = infer_gtr_impl(&counts, &InferGtrOptions::default())?;
  let alphabet = partition.read_arc().alphabet().clone();
  let W = Some(W);
  Ok(GTR::new(GTRParams { alphabet, mu, W, pi }))
}

#[derive(Clone, Debug)]
pub(super) struct MutationCounts {
  /// NxN matrix where each entry represents the observed number of transitions from state i to state j.
  pub nij: Array2<f64>,

  /// An N-vector representing the total time spent in each character state.
  pub Ti: Array1<f64>,

  /// An N-vector representing the state counts at the root of the phylogenetic tree.
  pub root_state: Array1<f64>,
}

#[derive(Clone, Debug, SmartDefault)]
pub(super) struct InferGtrOptions {
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

#[derive(Clone, Debug)]
pub(super) struct InferGtrResult {
  /// Substitution attempt matrix calculated during the GTR inference.
  pub W: Array2<f64>,
  /// Estimated equilibrium state frequencies.
  pub pi: Array1<f64>,
  /// Scale factor for the rates matrix.
  pub mu: f64,
}

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
pub(super) fn infer_gtr_impl(counts: &MutationCounts, options: &InferGtrOptions) -> Result<InferGtrResult, Report> {
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
    W.fill(0.0);
    W = (&(&nij.view() + &nij.t() + 2.0 * &pc_mat) / mu) / (&outer(&pi, Ti)? + &outer(Ti, &pi)? + 2.0 * &pc_mat);
    W /= avg_transition(&W, &pi)?;

    if fixed_pi.is_none() {
      pi = (&nij.sum_axis(Axis(1)) + &pc_mat.sum_axis(Axis(1)) + root_state)
        / (mu * W.dot(Ti) + root_state.sum() + pc_mat.sum_axis(Axis(1)));
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

pub(super) fn distance(pi_old: &Array1<f64>, pi: &Array1<f64>) -> f64 {
  (pi_old - pi).mapv(|x| x * x).sum().sqrt()
}

fn get_mutation_counts<P: PartitionWithGtrInference>(
  graph: &GraphAncestral,
  partition: &Arc<RwLock<P>>,
) -> Result<MutationCounts, Report> {
  let partition_guard = partition.read_arc();
  let alphabet = partition_guard.alphabet();

  let root_state = {
    let root = graph.get_exactly_one_root()?;
    let root_key = root.read_arc().key();
    let root_composition = partition_guard.get_seq_composition(root_key);
    Array1::<f64>::from_iter(
      alphabet
        .canonical()
        .map(|nuc| root_composition.get(nuc).unwrap_or(0) as f64),
    )
  };

  let N = alphabet.n_canonical();
  let mut nij = Array2::zeros((N, N));
  let mut Ti = Array1::zeros(N);

  for edge in graph.get_edges() {
    let edge_arc = edge.read_arc();
    let branch_length = edge_arc.payload().read_arc().branch_length().unwrap_or(0.0);
    let target_key = edge_arc.target();
    let edge_key = edge_arc.key();

    let node_composition = partition_guard.get_seq_composition(target_key);

    for (i, nuc) in alphabet.canonical().enumerate() {
      Ti[i] += branch_length * node_composition.get(nuc).unwrap_or(0) as f64;
    }

    let subs = partition_guard.get_edge_substitutions(edge_key, graph);
    for m in &subs {
      m.check_canonical(alphabet)?;
      let i = alphabet.index(m.qry());
      let j = alphabet.index(m.reff());
      nij[[i, j]] += 1.0;
      Ti[i] -= 0.5 * branch_length;
      Ti[j] += 0.5 * branch_length;
    }
  }

  Ok(MutationCounts { nij, Ti, root_state })
}
