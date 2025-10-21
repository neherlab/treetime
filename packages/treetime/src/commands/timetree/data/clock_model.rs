use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use eyre::Report;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ClockModel {
  pub mean_rate: Option<f64>,
  pub variance: Option<f64>,
  pub slack: Option<f64>,
}

pub fn infer_clock_model<N, E, D>(args: &TreetimeTimetreeArgs, graph: &Graph<N, E, D>) -> Result<ClockModel, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let _ = (args, graph);
  // TODO: implement initial clock rate estimation (e.g., root-to-tip regression) and capture uncertainty.
  todo!()
}

pub fn update_clock_model<N, E, D>(graph: &Graph<N, E, D>, model: &ClockModel) -> Result<ClockModel, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let _ = (graph, model);
  // TODO: iteratively refine the clock model using updated node time distributions.
  todo!()
}

pub fn clock_rate(model: &ClockModel) -> Option<f64> {
  model.mean_rate
}
