use crate::commands::timetree::coalescent::integration::compute_merger_rates;
use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_map::distribution_map;
use crate::distribution::distribution_multiplication::distribution_multiplication;
use crate::distribution::distribution_resample::distribution_resample;
use crate::graph::graph::GraphNodeForward;
use crate::graph::node::GraphNodeKey;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use eyre::{Context, Report};
use indexmap::IndexMap;
use std::sync::Arc;

/// Computes coalescent prior contributions for all nodes.
///
/// Returns map from node key to distribution that should be multiplied
/// with node's time distribution during backward pass.
///
/// Kingman coalescent probability density:
/// - For leaves: P(t) ∝ exp(-I(t)) where I(t) = ∫κ(t')dt' is cumulative merger rate
/// - For internal nodes with k children: P(t) ∝ λ(t)^(k-1) · exp(-I(t))
///   where λ(t) = k(k-1)/(2Tc) is total merger rate
pub fn compute_node_contributions(
  graph: &GraphAncestral,
  integral_merger_rate: &Distribution,
  tc_dist: &Distribution,
  lineage_counts: &Distribution,
  present_time: f64,
) -> Result<IndexMap<GraphNodeKey, Arc<Distribution>>, Report> {
  let mut contributions = IndexMap::new();

  graph.iter_breadth_first_forward(|node| {
    let contrib = compute_node_contribution_single(&node, integral_merger_rate, tc_dist, lineage_counts, present_time)
      .wrap_err_with(|| {
        format!(
          "When computing coalescent contributions for node (\"{}\") (#{})",
          node.payload.name.as_deref().unwrap_or(""),
          node.key,
        )
      })
      .unwrap();

    contributions.insert(node.key, Arc::new(contrib));
  });

  Ok(contributions)
}

fn compute_node_contribution_single(
  node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
  integral_merger_rate: &Distribution,
  tc_dist: &Distribution,
  lineage_counts: &Distribution,
  present_time: f64,
) -> Result<Distribution, Report> {
  let time_points = &node
    .payload
    .time_distribution
    .as_ref()
    .expect("Node missing time distribution")
    .t();

  let tbp_points = time_points.mapv(|t| present_time - t);
  let integral_at_node = distribution_resample(integral_merger_rate, &tbp_points)?;
  let exp_neg_integral = distribution_map(&integral_at_node, |x| (-x).exp())?;

  let contrib = if node.is_leaf {
    exp_neg_integral
  } else {
    let k_vals = lineage_counts.eval_many(&tbp_points)?;
    let tc_vals = tc_dist.eval_many(&tbp_points)?;

    let (_, total_rate) = compute_merger_rates(&k_vals, &tc_vals);
    let rate_dist = Distribution::function(tbp_points, total_rate)?;

    let n_children = node.child_edges.len();
    let multiplicity = (n_children - 1) as f64;
    let rate_power = distribution_map(&rate_dist, |r| r.powf(multiplicity))?;

    distribution_multiplication(&rate_power, &exp_neg_integral)?
  };

  Ok(contrib)
}
