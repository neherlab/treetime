use crate::coalescent::coalescent::CoalescentModel;
use crate::coalescent::edge_data::{coalescent_log_likelihood, collect_coalescent_edges};
use crate::coalescent::precomputed::CoalescentPrecomputed;
use crate::payload::traits::TimetreeNode;
use eyre::Report;
use treetime_distribution::Distribution;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

/// Computes the total coalescent log-likelihood of the tree under the given Tc.
///
/// For each non-root edge, evaluates the per-edge cost:
///   cost = I(t_merger) - I(t_node) - log(λ(t_merger)) * (m-1)/m
///
/// and returns LH = -Σ cost.
///
/// Multiplicity m is the number of children of the parent merger node. This corrects a v0 erratum where
/// `total_LH()` uses fixed multiplicity=2 for all edges (see
/// `docs/port-v0-errata/coalescent-total-lh-fixed-multiplicity.md`).
///
/// Accepts any `Distribution` for Tc (constant, skyline, or formula-based).
/// Nonconstant distributions are evaluated in decimal calendar years.
pub fn compute_coalescent_total_lh<N, E, D>(graph: &Graph<N, E, D>, tc_dist: &Distribution) -> Result<f64, Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let pre = CoalescentPrecomputed::from_graph(graph)?;
  let model = CoalescentModel::new(&pre, tc_dist)?;
  let edges = collect_coalescent_edges(graph)?;

  coalescent_log_likelihood(&edges, &model)
}
