use crate::commands::timetree::branch_distributions::BranchDistributions;
use crate::commands::timetree::date_constraints::DateConstraintSet;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use eyre::Report;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BackwardPassState {
  pub processed_edges: usize,
  pub iterations: usize,
}

pub fn run_backward_pass<N, E, D>(
  graph: &Graph<N, E, D>,
  constraints: &DateConstraintSet,
  distributions: &mut BranchDistributions,
) -> Result<BackwardPassState, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let _ = (graph, constraints, distributions);
  // TODO: perform post-order traversal and accumulate child-to-parent messages.
  todo!()
}
