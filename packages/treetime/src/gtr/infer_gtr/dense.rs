use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::infer_gtr::{InferGtrOptions, InferGtrResult, MutationCounts, infer_gtr_impl};
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::payload::ancestral::GraphAncestral;
use crate::representation::payload::dense::DenseSeqDis;
use eyre::Report;
use ndarray::{Array1, Array2};
use ndarray_stats::QuantileExt;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::HasBranchLength;
use treetime_primitives::Seq;

/// Infer GTR model from dense partition data.
pub fn infer_gtr_dense(partition: &Arc<RwLock<PartitionMarginalDense>>, graph: &GraphAncestral) -> Result<GTR, Report> {
  let counts = get_mutation_counts_dense(graph, partition)?;
  let InferGtrResult { W, pi, mu } = infer_gtr_impl(&counts, &InferGtrOptions::default())?;
  let alphabet = partition.read_arc().alphabet.clone();
  let W = Some(W);
  GTR::new(GTRParams { alphabet, mu, W, pi })
}

/// Get mutation counts from dense partition for GTR inference.
///
/// Extracts sequences from node profiles via argmax, then counts:
/// - `nij`: substitutions from state i to state j across all edges
/// - `Ti`: time spent in each state (branch length * composition)
/// - `root_state`: character composition at root
pub fn get_mutation_counts_dense(
  graph: &GraphAncestral,
  partition: &Arc<RwLock<PartitionMarginalDense>>,
) -> Result<MutationCounts, Report> {
  let partition = partition.read_arc();
  let alphabet = &partition.alphabet;
  let N = alphabet.n_canonical();

  let root_state = {
    let root = graph.get_exactly_one_root()?;
    let root_key = root.read_arc().key();
    let root_profile = &partition.nodes[&root_key].profile;
    let root_seq = prof2seq(root_profile, alphabet);
    count_composition(&root_seq, alphabet)
  };

  let mut nij = Array2::zeros((N, N));
  let mut Ti = Array1::zeros(N);

  for edge in graph.get_edges() {
    let edge_arc = edge.read_arc();
    let branch_length = edge_arc.payload().read_arc().branch_length().unwrap_or(0.0);
    let source_key = edge_arc.source();
    let target_key = edge_arc.target();

    let parent_profile = &partition.nodes[&source_key].profile;
    let child_profile = &partition.nodes[&target_key].profile;

    let parent_seq = prof2seq(parent_profile, alphabet);
    let child_seq = prof2seq(child_profile, alphabet);

    let child_composition = count_composition(&child_seq, alphabet);
    Ti += &(&child_composition * branch_length);

    count_substitutions(&parent_seq, &child_seq, alphabet, &mut nij, &mut Ti, branch_length)?;
  }

  Ok(MutationCounts { nij, Ti, root_state })
}

/// Extract sequence from profile by taking argmax at each position.
fn prof2seq(profile: &DenseSeqDis, alphabet: &Alphabet) -> Seq {
  let mut seq = Seq::new();
  for row in profile.dis.rows() {
    let argmax = row.argmax().unwrap_or(0);
    seq.push(alphabet.char(argmax));
  }
  seq
}

/// Count character composition of a sequence.
fn count_composition(seq: &Seq, alphabet: &Alphabet) -> Array1<f64> {
  let N = alphabet.n_canonical();
  let mut composition = Array1::zeros(N);
  for &ch in seq {
    if let Ok(idx) = alphabet.index(ch) {
      if idx < N {
        composition[idx] += 1.0;
      }
    }
  }
  composition
}

/// Count substitutions between parent and child sequences.
/// Also adjusts Ti for the branch length correction (same as sparse).
fn count_substitutions(
  parent_seq: &Seq,
  child_seq: &Seq,
  alphabet: &Alphabet,
  nij: &mut Array2<f64>,
  Ti: &mut Array1<f64>,
  branch_length: f64,
) -> Result<(), Report> {
  for (&parent_ch, &child_ch) in parent_seq.iter().zip(child_seq.iter()) {
    if parent_ch == child_ch {
      continue;
    }

    let parent_idx = alphabet.index(parent_ch);
    let child_idx = alphabet.index(child_ch);

    let (Ok(i), Ok(j)) = (parent_idx, child_idx) else {
      continue;
    };

    let N = alphabet.n_canonical();
    if i >= N || j >= N {
      continue;
    }

    nij[[j, i]] += 1.0;
    Ti[j] -= 0.5 * branch_length;
    Ti[i] += 0.5 * branch_length;
  }
  Ok(())
}
