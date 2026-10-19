use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::mean_field::NodeSurvival;
use crate::timetree::confidence::CI_FRACTION;
use crate::timetree::max_posterior_region::{MaxPosteriorRegion, max_posterior_region};
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct NodeBackward {
  pub subtree: TimeDistribution,
  pub to_parent: TimeDistribution,
}

pub(crate) fn log_marginal_likelihood(
  backward: &BTreeMap<GraphNodeKey, NodeBackward>,
  root: GraphNodeKey,
) -> Result<Option<f64>, Report> {
  backward
    .get(&root)
    .and_then(|node| node.subtree.as_deref())
    .map(Distribution::log_mass)
    .transpose()
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
  pub log_lh_time: Option<LogLh>,
  pub bad_branches: BTreeMap<GraphNodeKey, bool>,
  pub branches: BTreeMap<GraphEdgeKey, BranchLikelihood>,
  pub posterior: BTreeMap<GraphNodeKey, NodePosterior>,
  pub posterior_summaries: PosteriorSummaries,
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

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PosteriorSummaries {
  pub regions: bool,
  pub lineage_survival: bool,
}

pub(crate) type NodeTimes = BTreeMap<GraphNodeKey, Option<f64>>;

#[derive(Debug, Clone, PartialEq)]
pub struct BranchLikelihood {
  pub distribution: TimeDistribution,
  pub time_length: Option<f64>,
}

#[derive(Debug, Clone, Default, PartialEq)]
#[expect(
  clippy::partial_pub_fields,
  reason = "the survival function only feeds the coalescent prior inside this crate"
)]
pub struct NodePosterior {
  pub likely_time: Option<f64>,
  pub time: Option<f64>,
  pub contradicted: bool,
  pub support_bounds: Option<(f64, f64)>,
  pub region: Option<MaxPosteriorRegion>,
  pub(crate) survival: Option<Arc<NodeSurvival>>,
}

impl NodePosterior {
  pub(crate) fn summarize(
    key: GraphNodeKey,
    view: PosteriorView<'_>,
    summaries: PosteriorSummaries,
    bad_branch: bool,
  ) -> Result<Self, Report> {
    let distribution = view.distribution.map(Arc::as_ref);
    let region = match (distribution, view.time) {
      (Some(distribution), Some(_)) if summaries.regions => max_posterior_region(distribution, CI_FRACTION)
        .wrap_err_with(|| format!("When computing the confidence interval of node {key}"))?,
      _ => None,
    };
    let survival = match distribution {
      Some(distribution) if summaries.lineage_survival && !bad_branch => NodeSurvival::of_posterior(distribution)
        .wrap_err_with(|| format!("When computing the survival function of node {key}"))?
        .map(Arc::new),
      _ => None,
    };
    Ok(Self {
      likely_time: view.likely_time,
      time: view.time,
      contradicted: view.contradicted,
      support_bounds: distribution.and_then(Distribution::support_bounds),
      region,
      survival,
    })
  }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct PosteriorView<'a> {
  pub distribution: Option<&'a Arc<Distribution<NegLog>>>,
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
