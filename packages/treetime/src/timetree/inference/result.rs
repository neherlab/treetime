use crate::clock::date_constraints::DateConstraints;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;

#[derive(Debug, Clone, PartialEq)]
pub struct TimeBackward {
  pub subtree: BTreeMap<GraphNodeKey, TimeDistribution>,
  pub messages: BTreeMap<GraphEdgeKey, TimeDistribution>,
}

impl TimeBackward {
  /// Log marginal likelihood of the node times: the log mass of the root's subtree distribution, whose amplitude
  /// carries the levels of all messages and factors below it. Branch-length likelihoods enter peak-normalized, so the
  /// value omits their (rate- and root-independent) peaks. `None` when the root has no distribution.
  pub(crate) fn log_marginal_likelihood(&self, root: GraphNodeKey) -> Result<Option<f64>, Report> {
    self
      .subtree
      .get(&root)
      .and_then(Option::as_deref)
      .map(Distribution::log_mass)
      .transpose()
  }
}

pub(crate) fn given_times(graph: &Graph, constraints: &DateConstraints) -> Result<NodeTimes, Report> {
  graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      let time = likely_time_of_node(key, constraints.date_constraint(key).as_deref())?;
      Ok((key, time))
    })
    .collect()
}

pub(crate) fn likely_times(
  graph: &Graph,
  constraints: &DateConstraints,
  inference: &TimeInference,
) -> Result<NodeTimes, Report> {
  graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      let time = match constraints.date_constraint(key) {
        Some(constraint) => likely_time_of_node(key, Some(&constraint))?,
        None => inference.posterior[&key].likely_time,
      };
      Ok((key, time))
    })
    .collect()
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct TimeInference {
  /// Log marginal likelihood of the node times from the backward pass (see `fn TimeBackward::log_marginal_likelihood()`)
  pub log_lh_time: Option<LogLh>,
  pub bad_branches: BTreeMap<GraphNodeKey, bool>,
  pub branches: BTreeMap<GraphEdgeKey, BranchLikelihood>,
  pub posterior: BTreeMap<GraphNodeKey, NodePosterior>,
}

impl TimeInference {
  #[must_use]
  pub(crate) fn node_times(&self) -> NodeTimes {
    self
      .posterior
      .iter()
      .map(|(key, posterior)| (*key, posterior.time))
      .collect()
  }
}

pub(crate) type NodeTimes = BTreeMap<GraphNodeKey, Option<f64>>;

#[derive(Debug, Clone, PartialEq)]
pub struct BranchLikelihood {
  pub distribution: TimeDistribution,
  pub time_length: Option<f64>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct NodePosterior {
  pub distribution: TimeDistribution,
  /// Point estimate of the node's time: the posterior peak, and for the root the posterior median
  pub likely_time: Option<f64>,
  pub time: Option<f64>,
  pub contradicted: bool,
}

pub type TimeDistribution = Option<Arc<Distribution<NegLog>>>;

pub(crate) fn likely_time_of_node(
  key: GraphNodeKey,
  distribution: Option<&Distribution<NegLog>>,
) -> Result<Option<f64>, Report> {
  distribution.map_or(Ok(None), |distribution| {
    distribution
      .likely_time()
      .wrap_err_with(|| format!("When finding the most likely time of node {key}"))
  })
}
