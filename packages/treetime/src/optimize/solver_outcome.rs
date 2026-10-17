use crate::node_label::node_label;
use crate::optimize::params::BranchOptMethod;
use crate::progress::{LogSink, RunWarning, RunWarningKind};
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct EdgeSolution {
  pub(crate) branch_length: f64,
  pub(crate) convergence: Convergence,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Convergence {
  Converged,
  IterationLimit,
  SearchBound,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct UnconvergedEdges {
  pub(crate) iteration_limit: Vec<GraphEdgeKey>,
  pub(crate) search_bound: Vec<GraphEdgeKey>,
}

impl UnconvergedEdges {
  pub(crate) fn of(solutions: &[(GraphEdgeKey, EdgeSolution)]) -> Self {
    let edges_with = |convergence: Convergence| {
      solutions
        .iter()
        .filter(|(_, solution)| solution.convergence == convergence)
        .map(|&(edge_key, _)| edge_key)
        .collect_vec()
    };
    Self {
      iteration_limit: edges_with(Convergence::IterationLimit),
      search_bound: edges_with(Convergence::SearchBound),
    }
  }

  fn is_empty(&self) -> bool {
    self.iteration_limit.is_empty() && self.search_bound.is_empty()
  }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct ConvergenceSummary {
  rounds: usize,
  unconverged_rounds: usize,
  iteration_limit: BTreeSet<String>,
  search_bound: BTreeSet<String>,
}

impl ConvergenceSummary {
  pub(crate) fn record(
    &mut self,
    round: &UnconvergedEdges,
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
  ) {
    self.rounds += 1;
    if round.is_empty() {
      return;
    }
    self.unconverged_rounds += 1;
    let label = |edge_key: &GraphEdgeKey| node_label(names, graph[*edge_key].target());
    self.iteration_limit.extend(round.iteration_limit.iter().map(label));
    self.search_bound.extend(round.search_bound.iter().map(label));
  }

  pub(crate) fn merge(&mut self, other: Self) {
    self.rounds += other.rounds;
    self.unconverged_rounds += other.unconverged_rounds;
    self.iteration_limit.extend(other.iteration_limit);
    self.search_bound.extend(other.search_bound);
  }

  pub(crate) fn report(&self, method: BranchOptMethod, log: &dyn LogSink) {
    if self.unconverged_rounds == 0 {
      return;
    }
    let names = self.iteration_limit.union(&self.search_bound).cloned().collect_vec();
    let causes = [
      (self.iteration_limit.len(), "it reached its iteration limit on"),
      (
        self.search_bound.len(),
        "it stopped at the upper end of its search interval, where the likelihood still increased, on",
      ),
    ]
    .into_iter()
    .filter(|&(count, _)| count > 0)
    .map(|(count, cause)| format!("{cause} {}", counted(count, "edge")))
    .join("; ");
    log.warning(&RunWarning {
      kind: RunWarningKind::BranchLengthsNotConverged,
      message: format!(
        "The {method} branch-length optimizer did not converge on {} in {} of {}: {causes}. These edges keep the best length the optimizer evaluated, which may not maximize the likelihood",
        counted(names.len(), "edge"),
        self.unconverged_rounds,
        counted(self.rounds, "optimization round"),
      ),
      names,
    });
  }
}

fn counted(count: usize, noun: &str) -> String {
  if count == 1 {
    format!("{count} {noun}")
  } else {
    format!("{count} {noun}s")
  }
}
