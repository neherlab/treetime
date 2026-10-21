use crate::evolve::validate_edge_duration;
use crate::substitution_model::SubstitutionModel;
use eyre::Report;
use itertools::izip;
use rand::Rng;
use rand_distr::{Distribution, LogNormal};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_error, make_report};

pub fn exact(
  durations: &BTreeMap<GraphEdgeKey, f64>,
  model: &SubstitutionModel,
) -> Result<BTreeMap<GraphEdgeKey, f64>, Report> {
  durations
    .iter()
    .map(|(&edge, &duration)| {
      validate_edge_duration(edge, duration)?;
      Ok((edge, model.mu() * duration))
    })
    .collect()
}

pub fn realized(graph: &Graph, sequences: &BTreeMap<GraphNodeKey, Seq>) -> Result<BTreeMap<GraphEdgeKey, f64>, Report> {
  graph
    .get_edges()
    .map(|edge| {
      let parent = &sequences[&edge.source()];
      let child = &sequences[&edge.target()];
      if parent.len() != child.len() {
        return make_error!(
          "Sequences of nodes {} and {} have different lengths: {} and {}",
          edge.source(),
          edge.target(),
          parent.len(),
          child.len()
        );
      }
      if parent.is_empty() {
        return make_error!("Sequence of node {} is empty", edge.source());
      }
      let differences = izip!(parent.iter(), child.iter())
        .filter(|(left, right)| left != right)
        .count();
      Ok((edge.key(), differences.to_f64_exact() / parent.len().to_f64_exact()))
    })
    .collect()
}

pub fn perturbed(
  lengths: &BTreeMap<GraphEdgeKey, f64>,
  sigma: f64,
  rng: &mut impl Rng,
) -> Result<BTreeMap<GraphEdgeKey, f64>, Report> {
  if !(sigma.is_finite() && sigma >= 0.0) {
    return make_error!("Perturbation sigma must be finite and non-negative, got {sigma}");
  }
  let noise =
    LogNormal::new(0.0, sigma).map_err(|error| make_report!("Log-normal noise with sigma {sigma}: {error}"))?;
  Ok(
    lengths
      .iter()
      .map(|(&edge, &length)| (edge, length * noise.sample(rng)))
      .collect(),
  )
}
