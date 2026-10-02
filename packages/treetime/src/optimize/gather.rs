use crate::partition::marginal::reconstruction::MarginalReconstruction;
use crate::partition::optimize::contribution::OptimizationContribution;
use eyre::Report;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

pub(crate) fn gather_edge_contributions(
  graph: &Graph,
  reconstruction: &MarginalReconstruction,
) -> Result<BTreeMap<GraphEdgeKey, Vec<OptimizationContribution>>, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(
      |edge_ref| -> Result<(GraphEdgeKey, Vec<OptimizationContribution>), Report> {
        let edge_key = edge_ref.key();
        Ok((edge_key, vec![reconstruction.create_edge_contribution(edge_key)?]))
      },
    )
    .collect()
}

pub(crate) fn gather_edge_indel_counts(
  graph: &Graph,
  reconstruction: &MarginalReconstruction,
) -> BTreeMap<GraphEdgeKey, usize> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| {
      let edge_key = edge_ref.key();
      (edge_key, reconstruction.edge_indel_count(edge_key))
    })
    .collect()
}

pub(crate) fn gather_edge_sub_counts(
  graph: &Graph,
  reconstruction: &MarginalReconstruction,
) -> Result<BTreeMap<GraphEdgeKey, usize>, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| -> Result<(GraphEdgeKey, usize), Report> {
      let edge_key = edge_ref.key();
      Ok((edge_key, reconstruction.edge_subs(graph, edge_key)?.len()))
    })
    .collect()
}

pub(crate) fn gather_edge_effective_lengths(
  graph: &Graph,
  reconstruction: &MarginalReconstruction,
) -> Result<BTreeMap<GraphEdgeKey, usize>, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| -> Result<(GraphEdgeKey, usize), Report> {
      let edge_key = edge_ref.key();
      Ok((edge_key, reconstruction.edge_effective_length(graph, edge_key)?))
    })
    .collect()
}
