use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::{
  InferGtrOptions, InferGtrResult, MutationCounts, count_branch_mutations, infer_gtr_impl,
};
use crate::partition::fitch::partition::PartitionFitch;
use crate::progress::LogSink;
use crate::seq::mutation::Sub;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

pub(crate) fn infer_gtr_fitch(
  partition: &PartitionFitch,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<GTR, Report> {
  let counts = get_mutation_counts_fitch(graph, partition, branch_lengths)?;
  let InferGtrResult { W, pi, mu } = infer_gtr_impl(&counts, &InferGtrOptions::default(), log)?;
  let n_states = partition.alphabet.n_canonical();
  let W = Some(W);
  GTR::builder().n_states(n_states).mu(mu).maybe_W(W).pi(pi).build()
}

pub(crate) fn get_mutation_counts_fitch(
  graph: &Graph,
  partition: &PartitionFitch,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<MutationCounts, Report> {
  let root_key = graph.get_exactly_one_root()?.key();
  let empty: &[Sub] = &[];
  let branches = graph.get_edges().map(|edge| {
    let edge_key = edge.key();
    let subs = partition.edges.get(&edge_key).map_or(empty, |e| e.fitch_subs());
    (
      branch_lengths[&edge_key],
      &partition.nodes[&edge.target()].seq.composition,
      subs,
    )
  });
  count_branch_mutations(
    &partition.alphabet,
    &partition.nodes[&root_key].seq.composition,
    branches,
  )
}
