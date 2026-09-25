use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::node_time::{CoalescentNodeTime, CoalescentNodeTimes};
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

#[derive(Debug, Clone, PartialEq)]
pub struct TimeInference {
  pub bad_branches: BTreeMap<GraphNodeKey, bool>,
  pub branches: BTreeMap<GraphEdgeKey, BranchLikelihood>,
  pub backward: TimeBackward,
  pub posterior: BTreeMap<GraphNodeKey, NodePosterior>,
}

impl TimeInference {
  #[must_use]
  pub(crate) fn node_times(&self) -> BTreeMap<GraphNodeKey, Option<f64>> {
    self
      .posterior
      .iter()
      .map(|(key, posterior)| (*key, posterior.time))
      .collect()
  }

  pub(crate) fn coalescent_node_times(&self) -> Result<CoalescentNodeTimes, Report> {
    self
      .posterior
      .iter()
      .map(|(key, posterior)| {
        let entry = CoalescentNodeTime {
          time: posterior.time,
          time_dist_likely: likely_time_of_node(*key, posterior.distribution.as_deref())?,
          bad_branch: self.bad_branches[key],
        };
        Ok((*key, entry))
      })
      .collect()
  }
}

pub type TimeMessage = Option<Arc<Distribution<NegLog>>>;

#[derive(Debug, Clone, PartialEq)]
pub struct BranchLikelihood {
  pub distribution: Option<Arc<Distribution<NegLog>>>,
  pub time_length: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct TimeBackward {
  pub subtree: BTreeMap<GraphNodeKey, TimeMessage>,
  pub messages: BTreeMap<GraphEdgeKey, TimeMessage>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct NodePosterior {
  pub distribution: Option<Arc<Distribution<NegLog>>>,
  pub time: Option<f64>,
  pub contradicted: bool,
}

pub(crate) fn likely_times(
  graph: &Graph,
  constraints: &DateConstraints,
  inference: Option<&TimeInference>,
) -> Result<BTreeMap<GraphNodeKey, Option<f64>>, Report> {
  graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      let distribution = constraints
        .date_constraint(key)
        .or_else(|| inference.and_then(|inference| inference.posterior[&key].distribution.clone()));
      let time = likely_time_of_node(key, distribution.as_deref())?;
      Ok((key, time))
    })
    .collect()
}

#[must_use]
pub(crate) fn unit_gammas(graph: &Graph) -> BTreeMap<GraphEdgeKey, f64> {
  graph.get_edges().map(|edge| (edge.key(), 1.0)).collect()
}

fn likely_time_of_node(key: GraphNodeKey, distribution: Option<&Distribution<NegLog>>) -> Result<Option<f64>, Report> {
  distribution.map_or(Ok(None), |distribution| {
    distribution
      .likely_time()
      .wrap_err_with(|| format!("When finding the most likely time of node {key}"))
  })
}
