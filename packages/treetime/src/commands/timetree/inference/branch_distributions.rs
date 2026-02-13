use crate::distribution::distribution::Distribution;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub enum BranchDistributionContext {
  InputPrior,
  Marginal,
  Joint,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BranchDistribution {
  pub context: BranchDistributionContext,
  pub profile: Arc<Distribution>, // TODO: revisit storage if mutable distributions become necessary.
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct BranchDistributions {
  pub per_edge: BTreeMap<GraphEdgeKey, BranchDistribution>,
}

pub fn build_branch_distributions<N, E, D>(
  graph: &Graph<N, E, D>,
  context: BranchDistributionContext,
) -> Result<BranchDistributions, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let _ = (graph, context);
  // TODO: derive edge-wise time priors from branch lengths or marginal reconstruction messages.
  todo!()
}

pub fn branch_distribution(distributions: &BranchDistributions, key: GraphEdgeKey) -> Option<&BranchDistribution> {
  distributions.per_edge.get(&key)
}
