use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, MutationCounts, infer_gtr_impl};
use crate::partition::fitch::partition::PartitionFitch;
use crate::seq::mutation::Sub;
use eyre::Report;
use ndarray::{Array1, Array2};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

pub(crate) fn infer_gtr_fitch(
  partition: &PartitionFitch,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<GTR, Report> {
  let counts = get_mutation_counts_fitch(graph, partition, branch_lengths)?;
  let InferGtrResult { W, pi, mu } = infer_gtr_impl(&counts, &InferGtrOptions::default())?;
  let n_states = partition.alphabet.n_canonical();
  let W = Some(W);
  GTR::builder().n_states(n_states).mu(mu).maybe_W(W).pi(pi).build()
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn get_mutation_counts_fitch(
  graph: &Graph,
  partition: &PartitionFitch,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<MutationCounts, Report> {
  let alphabet = &partition.alphabet;

  let root_state = {
    let root = graph.get_exactly_one_root()?;
    let root_key = root.key();
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
    let edge_arc = edge;
    let target_key = edge_arc.target();
    let edge_key = edge_arc.key();
    let branch_length = branch_lengths[&edge_key];

    let node_composition = &partition.nodes[&target_key].seq.composition;

    for (i, nuc) in alphabet.canonical().enumerate() {
      Ti[i] += branch_length * node_composition.get(nuc).unwrap_or(0) as f64;
    }

    let empty: &[Sub] = &[];
    let subs = partition.edges.get(&edge_key).map_or(empty, |e| e.fitch_subs());
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
