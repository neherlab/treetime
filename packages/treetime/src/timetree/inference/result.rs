use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::node_time::{CoalescentNodeTime, CoalescentNodeTimes};
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct TimeInference {
  pub bad_branches: BTreeMap<GraphNodeKey, bool>,
  pub branches: BTreeMap<GraphEdgeKey, BranchLikelihood>,
  pub posterior: BTreeMap<GraphNodeKey, NodePosterior>,
}

pub(crate) type NodeTimes = BTreeMap<GraphNodeKey, Option<f64>>;

impl TimeInference {
  #[must_use]
  pub(crate) fn node_times(&self) -> NodeTimes {
    self
      .posterior
      .iter()
      .map(|(key, posterior)| (*key, posterior.time))
      .collect()
  }

  pub(crate) fn coalescent_node_times(&self) -> CoalescentNodeTimes {
    self
      .posterior
      .iter()
      .map(|(key, posterior)| {
        let entry = CoalescentNodeTime {
          time: posterior.time,
          time_dist_likely: posterior.likely_time,
          bad_branch: self.bad_branches[key],
        };
        (*key, entry)
      })
      .collect()
  }
}

pub type TimeDistribution = Option<Arc<Distribution<NegLog>>>;

#[derive(Debug, Clone, PartialEq)]
pub struct BranchLikelihood {
  pub distribution: TimeDistribution,
  pub time_length: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct TimeBackward {
  pub subtree: BTreeMap<GraphNodeKey, TimeDistribution>,
  pub messages: BTreeMap<GraphEdgeKey, TimeDistribution>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct NodePosterior {
  pub distribution: TimeDistribution,
  pub likely_time: Option<f64>,
  pub time: Option<f64>,
  pub contradicted: bool,
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
