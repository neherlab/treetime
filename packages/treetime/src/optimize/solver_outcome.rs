use crate::node_label::node_label;
use crate::optimize::params::BranchOptMethod;
use crate::progress::{LogSink, RunWarning, RunWarningKind};
use itertools::Itertools;
use std::collections::BTreeMap;
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
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct UnconvergedEdges {
  pub(crate) edges: Vec<GraphEdgeKey>,
  pub(crate) optimized: usize,
}

impl UnconvergedEdges {
  pub(crate) fn report(
    &self,
    method: BranchOptMethod,
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    log: &dyn LogSink,
  ) {
    if self.edges.is_empty() {
      return;
    }
    let names = self
      .edges
      .iter()
      .map(|&edge_key| node_label(names, graph[edge_key].target()))
      .sorted()
      .dedup()
      .collect_vec();
    log.warning(&RunWarning {
      kind: RunWarningKind::BranchLengthsNotConverged,
      message: format!(
        "The {method} branch-length optimizer reached its iteration limit on {} of {} edges. These edges keep the best length the optimizer evaluated, which may not maximize the likelihood",
        self.edges.len(),
        self.optimized
      ),
      names,
    });
  }
}
