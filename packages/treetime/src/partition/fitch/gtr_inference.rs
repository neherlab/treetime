use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, MutationCounts, infer_gtr_from_counts};
use crate::partition::fitch::partition::PartitionFitch;
use crate::progress::LogSink;
use crate::seq::mutation::Sub;
use eyre::Report;
use ndarray::{Array1, Array2};
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::schedule::GraphSchedule;
use treetime_utils::cast::ToF64Exact;

pub(crate) fn infer_gtr_fitch(
  partition: &PartitionFitch,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  log: &dyn LogSink,
) -> Result<GTR, Report> {
  let counts = get_mutation_counts_fitch(graph, partition, branch_lengths)?;
  let InferGtrResult { W, pi, mu } = infer_gtr_from_counts(&counts, &InferGtrOptions::default(), log)?;
  let n_states = partition.alphabet.n_canonical();
  let W = Some(W);
  GTR::builder().n_states(n_states).mu(mu).maybe_W(W).pi(pi).build()
}

pub(crate) fn get_mutation_counts_fitch(
  graph: &Graph,
  partition: &PartitionFitch,
  branch_lengths: &EdgeMap<f64>,
) -> Result<MutationCounts, Report> {
  let alphabet = &partition.alphabet;

  let schedule = GraphSchedule::new(graph)?;
  let mut root_state = Array1::<f64>::zeros(alphabet.n_canonical());
  for root_key in partition.inheritance.roots(&schedule) {
    let root_composition = &partition.nodes[&root_key].seq.composition;
    for (count, nuc) in root_state.iter_mut().zip(alphabet.canonical()) {
      *count += root_composition.get(nuc).unwrap_or(0).to_f64_exact();
    }
  }

  let N = alphabet.n_canonical();
  let mut nij = Array2::zeros((N, N));
  let mut Ti = Array1::zeros(N);

  for edge in graph.get_edges() {
    let target_key = edge.target();
    let edge_key = edge.key();
    if !partition.inheritance.carries(edge_key) {
      continue;
    }
    let branch_length = branch_lengths[&edge_key];

    let node_composition = &partition.nodes[&target_key].seq.composition;

    for (i, nuc) in alphabet.canonical().enumerate() {
      Ti[i] += branch_length * node_composition.get(nuc).unwrap_or(0).to_f64_exact();
    }

    let empty: &[Sub] = &[];
    let subs = partition.edges.get(edge_key).map_or(empty, |e| e.fitch_subs());
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
