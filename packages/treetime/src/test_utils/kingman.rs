use eyre::Report;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_simulation::coalescent::SimulatedTree;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::make_error;

pub(crate) fn kingman_neg_log_density(nodes: &[NodeAge], tc: f64) -> Result<f64, Report> {
  if !(tc.is_finite() && tc > 0.0) {
    return make_error!("Kingman density needs a positive finite Tc, got {tc}");
  }
  let exposure = pair_exposure(nodes)?;
  Ok(exposure.mergers.to_f64_exact() * tc.ln() + exposure.pair_time / tc)
}

pub(crate) fn pair_exposure(nodes: &[NodeAge]) -> Result<PairExposure, Report> {
  let exposure = nodes
    .iter()
    .sorted_by_key(|node| (OrderedFloat(node.age), node.n_children))
    .fold(PairExposure::default(), PairExposure::advance);
  if exposure.lineages != 1 {
    return make_error!(
      "Node ages leave {} lineages after the oldest node; a tree ends with one",
      exposure.lineages
    );
  }
  Ok(exposure)
}

pub(crate) fn node_ages_from_dates(graph: &Graph, dates: &BTreeMap<GraphNodeKey, f64>) -> Vec<NodeAge> {
  graph
    .get_nodes()
    .map(|node| NodeAge {
      age: -dates[&node.key()],
      n_children: node.outbound().len(),
    })
    .collect()
}

pub(crate) fn simulated_node_ages(tree: &SimulatedTree) -> Vec<NodeAge> {
  tree
    .graph
    .get_nodes()
    .map(|node| NodeAge {
      age: tree.tau[&node.key()],
      n_children: node.outbound().len(),
    })
    .collect()
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct NodeAge {
  pub age: f64,
  pub n_children: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct PairExposure {
  pub pair_time: f64,
  pub mergers: usize,
  age: f64,
  lineages: isize,
}

impl PairExposure {
  fn advance(self, node: &NodeAge) -> Self {
    let lineages = self.lineages.to_f64_exact();
    let pair_time = self.pair_time + lineages * (lineages - 1.0) / 2.0 * (node.age - self.age);
    Self {
      pair_time,
      mergers: self.mergers + node.n_children.saturating_sub(1),
      age: node.age,
      lineages: self.lineages + 1 - isize::try_from(node.n_children).expect("child count must fit into isize"),
    }
  }
}
