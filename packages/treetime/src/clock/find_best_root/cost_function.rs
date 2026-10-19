use crate::branch_lengths::branch_length_or_zero;
use crate::clock::clock_set::ClockSet;
use crate::clock::clock_state::ClockState;
use crate::clock::find_best_root::params::RootObjective;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;

pub(crate) struct BranchPointCostFunction {
  to_parent: ClockSet,
  to_child: ClockSet,
  branch_length: f64,
  objective: RootObjective,
}

impl BranchPointCostFunction {
  pub(crate) fn new(
    state: &ClockState,
    edge: GraphEdgeKey,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    objective: RootObjective,
  ) -> Self {
    let edge_state = state.edge(edge);
    Self {
      to_parent: edge_state.clock_to_parent.clone(),
      to_child: edge_state.clock_to_child.clone(),
      branch_length: branch_length_or_zero(branch_lengths, edge),
      objective,
    }
  }

  pub(crate) fn evaluate_clock_set(&self, x: f64) -> ClockSet {
    self.to_child.shift(self.branch_length * x) + self.to_parent.shift(self.branch_length * (1.0 - x))
  }

  pub(crate) fn score_clock_set(&self, clock_set: &ClockSet) -> f64 {
    self.objective.score(clock_set)
  }
}
