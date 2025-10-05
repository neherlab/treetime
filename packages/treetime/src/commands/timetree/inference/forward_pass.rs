use crate::commands::timetree::data::date_constraints::DateConstraintSet;
use crate::commands::timetree::inference::backward_pass::BackwardPassState;
use crate::commands::timetree::inference::branch_distributions::BranchDistributions;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use eyre::Report;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ForwardPassState {
  pub processed_edges: usize,
  pub iterations: usize,
}

pub fn run_forward_pass<N, E, D>(
  graph: &Graph<N, E, D>,
  constraints: &DateConstraintSet,
  distributions: &mut BranchDistributions,
  backward: &BackwardPassState,
) -> Result<ForwardPassState, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let _ = (graph, constraints, distributions, backward);
  // TODO: run preorder traversal to propagate parent messages and form final posteriors.
  todo!()
}
