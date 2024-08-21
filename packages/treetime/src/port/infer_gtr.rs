use crate::alphabet::alphabet::Alphabet;
use crate::graph::edge::Weighted;
use crate::gtr::gtr::avg_transition;
use crate::port::seq_sparse::SparseGraph;
use crate::utils::ndarray::outer;
use eyre::Report;
use log::warn;
use ndarray::{Array1, Array2, Axis};
use smart_default::SmartDefault;

#[derive(Clone, Debug)]
pub struct MutationCounts {
  /// NxN matrix where each entry represents the observed number of transitions from state i to state j.
  pub nij: Array2<f64>,

  /// An N-vector representing the total time spent in each character state.
  pub Ti: Array1<f64>,

  /// An N-vector representing the state counts at the root of the phylogenetic tree.
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

#[derive(Clone, Debug)]
pub struct InferGtrResult {
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
pub fn infer_gtr(counts: &MutationCounts, options: &InferGtrOptions) -> Result<InferGtrResult, Report> {
  let MutationCounts { nij, Ti, root_state } = counts;
  let InferGtrOptions {
    fixed_pi,
    pc,
    dp,
    max_iter,
  } = options;

  let N = Ti.len();
  let pc_mat = Array2::from_elem((N, N), *pc);
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

fn distance(pi_old: &Array1<f64>, pi: &Array1<f64>) -> f64 {
  (pi_old - pi).mapv(|x| x * x).sum().sqrt()
}

pub fn get_mutation_counts(graph: &SparseGraph, alphabet: &Alphabet) -> Result<MutationCounts, Report> {
  let root = graph.get_exactly_one_root()?.read_arc().payload().read_arc();

  let seq = &root.sparse_partitions[0].seq;
  let root_state: Array1<f64> = Array1::from_iter(
    alphabet
      .canonical()
      .map(|nuc| seq.composition.get(nuc).unwrap_or(0) as f64),
  );

  let N = alphabet.n_canonical();
  let mut nij = Array2::zeros((N, N));
  let mut Ti = Array1::zeros(N);
  for edge in graph.get_edges() {
    let target_seq = &graph
      .get_node(edge.read_arc().target())
      .unwrap()
      .read_arc()
      .payload()
      .read_arc()
      .sparse_partitions[0]
      .seq;

    let edge = edge.read_arc().payload().read_arc();
    let branch_length = edge.weight().unwrap_or(0.0);

    for (i, nuc) in alphabet.canonical().enumerate() {
      Ti[i] += branch_length * target_seq.composition.get(nuc).unwrap_or(0) as f64;
    }

    for m in &edge.sparse_partitions[0].muts {
      let i = alphabet.index(m.qry);
      let j = alphabet.index(m.reff);
      nij[[i, j]] += 1.0;
      Ti[i] -= 0.5 * branch_length;
      Ti[j] += 0.5 * branch_length;
    }
  }

  Ok(MutationCounts { nij, Ti, root_state })
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::alphabet::alphabet::AlphabetName;
  use crate::gtr::gtr::{GTRParams, GTR};
  use crate::io::fasta::read_many_fasta_str;
  use crate::io::nwk::nwk_read_str;
  use crate::port::fitch::compress_sequences;
  use crate::port::seq_partitions::SeqPartition;
  use crate::pretty_assert_ulps_eq;
  use indoc::indoc;
  use ndarray::array;

  #[test]
  fn test_infer_gtr_only() -> Result<(), Report> {
    let nij = array![
      [0.0, 1.0, 2.0, 1.0],
      [1.0, 0.0, 3.0, 2.0],
      [2.0, 3.0, 0.0, 1.0],
      [2.0, 3.0, 3.0, 0.0]
    ];
    let Ti = array![12.0, 20.0, 14.0, 12.4];
    let root_state = array![3.0, 2.0, 3.0, 4.0];

    let actual = infer_gtr(
      &MutationCounts { nij, Ti, root_state },
      &InferGtrOptions {
        pc: 0.1,
        ..InferGtrOptions::default()
      },
    )?;

    pretty_assert_ulps_eq!(
      array![0.20908015, 0.24528811, 0.20925859, 0.33637315],
      &actual.pi,
      epsilon = 1e-3
    );

    pretty_assert_ulps_eq!(0.4004706866848001, actual.mu, epsilon = 1e-3);

    pretty_assert_ulps_eq!(1.0, avg_transition(&actual.W, &actual.pi)?, epsilon = 1e-5);

    Ok(())
  }

  #[test]
  fn test_infer_gtr_with_mutation_counts() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;
    let inputs = read_many_fasta_str(indoc! {r#"
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
      "#})?;

    let graph: SparseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
    let gtr = GTR::new(GTRParams {
      alphabet: alphabet.clone(),
      mu: 1.0,
      W: None,
      pi: array![0.2, 0.3, 0.15, 0.35],
    })?;

    let partitions = vec![SeqPartition::new(gtr, inputs)?];
    compress_sequences(&graph, &partitions)?;

    let counts_actual = get_mutation_counts(&graph, &alphabet)?;
    let counts_expected = MutationCounts {
      nij: array![[0., 0., 0., 0.], [2., 0., 0., 1.], [3., 2., 0., 0.], [0., 1., 1., 0.]],
      Ti: array![1.98, 2.945, 2.515, 2.64],
      root_state: array![4.0, 3.0, 3.0, 4.0],
    };

    pretty_assert_ulps_eq!(counts_expected.nij, counts_actual.nij, epsilon = 1e-9);
    pretty_assert_ulps_eq!(counts_expected.Ti, counts_actual.Ti, epsilon = 1e-9);
    pretty_assert_ulps_eq!(counts_expected.root_state, counts_actual.root_state, epsilon = 1e-9);

    let actual = infer_gtr(
      &counts_actual,
      &InferGtrOptions {
        pc: 0.1,
        ..InferGtrOptions::default()
      },
    )?;

    let expected = InferGtrResult {
      W: array![
        [0.0, 2.1751124, 2.95601658, 0.18620301],
        [2.1751124, 0.0, 1.40528091, 1.41465696],
        [2.95601658, 1.40528091, 0.0, 0.74490315],
        [0.18620301, 1.41465696, 0.74490315, 0.0]
      ],
      pi: array![0.14878846, 0.24051536, 0.31239203, 0.29830414],
      mu: 0.9471364432348814,
    };

    pretty_assert_ulps_eq!(expected.W, actual.W, epsilon = 1e-9);
    pretty_assert_ulps_eq!(expected.pi, actual.pi, epsilon = 1e-9);
    pretty_assert_ulps_eq!(expected.mu, actual.mu, epsilon = 1e-9);

    Ok(())
  }
}
