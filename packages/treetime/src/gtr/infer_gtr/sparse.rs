use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, MutationCounts, infer_gtr_impl};
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use ndarray::{Array1, Array2};
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::HasBranchLength;

/// Infer GTR model from sparse partition data.
pub fn infer_gtr_sparse(
  partition: &Arc<RwLock<PartitionMarginalSparse>>,
  graph: &GraphAncestral,
) -> Result<GTR, Report> {
  let counts = get_mutation_counts_sparse(graph, partition)?;
  let InferGtrResult { W, pi, mu } = infer_gtr_impl(&counts, &InferGtrOptions::default())?;
  let n_states = partition.read_arc().alphabet.n_canonical();
  let W = Some(W);
  GTR::new(GTRParams { n_states, mu, W, pi })
}

pub fn get_mutation_counts_sparse(
  graph: &GraphAncestral,
  partition: &Arc<RwLock<PartitionMarginalSparse>>,
) -> Result<MutationCounts, Report> {
  let partition = partition.read_arc();
  let alphabet = &partition.alphabet;

  let root_state = {
    let root = graph.get_exactly_one_root()?;
    let root_key = root.read_arc().key();
    let root_composition = &partition.nodes[&root_key].seq.composition;
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

    let node_composition = &partition.nodes[&target_key].seq.composition;

    for (i, nuc) in alphabet.canonical().enumerate() {
      Ti[i] += branch_length * node_composition.get(nuc).unwrap_or(0) as f64;
    }

    let subs = &partition.edges[&edge_key].subs;
    for m in subs {
      m.check_canonical(alphabet)?;
      let i = alphabet.index(m.qry())?;
      let j = alphabet.index(m.reff())?;
      nij[[i, j]] += 1.0;
      Ti[i] -= 0.5 * branch_length;
      Ti[j] += 0.5 * branch_length;
    }
  }

  Ok(MutationCounts { nij, Ti, root_state })
}
