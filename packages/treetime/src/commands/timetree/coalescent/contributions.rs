use crate::commands::timetree::coalescent::integration::compute_merger_rates;
use crate::commands::timetree::coalescent::piecewise_constant_fn::PiecewiseConstantFn;
use crate::commands::timetree::coalescent::piecewise_linear_fn::PiecewiseLinearFn;
use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::graph::graph::GraphNodeForward;
use crate::graph::node::GraphNodeKey;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use eyre::{Context, Report};
use indexmap::IndexMap;
use ndarray::Array1;
use std::sync::Arc;

/// Computes coalescent prior contributions for all nodes.
///
/// Returns map from node key to distribution that should be multiplied
/// with node's time distribution during backward pass.
///
/// Kingman coalescent contribution formula:
/// - For leaves: exp(I(t)) where I(t) = ∫₀ᵗ κ(t')dt' is cumulative merger rate
///   (stored as -I(t) in neg-log space, which converts to exp(I(t)) in probability)
/// - For internal nodes with k children: λ(t)^m · exp(-m·I(t))
///   where λ(t) = k(t)·(k(t)-1)/(2·Tc(t)) is total merger rate
///   and m = k - 1 is multiplicity (number of mergers)
pub fn compute_node_contributions(
  graph: &GraphAncestral,
  integral_merger_rate: &PiecewiseLinearFn,
  tc_dist: &Distribution,
  lineage_counts: &PiecewiseConstantFn,
  present_time: f64,
) -> Result<IndexMap<GraphNodeKey, Arc<Distribution>>, Report> {
  let mut contributions = IndexMap::new();

  graph.iter_breadth_first_forward(|node| {
    // Skip leaf nodes - they represent observed samples at known times,
    // not coalescence events, so they should not have coalescent priors
    if node.is_leaf {
      return;
    }

    // Skip nodes without time distributions
    if node.payload.time_distribution.is_none() {
      return;
    }

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
  integral_merger_rate: &PiecewiseLinearFn,
  tc_dist: &Distribution,
  lineage_counts: &PiecewiseConstantFn,
  _present_time: f64,
) -> Result<Distribution, Report> {
  // Use the integral_merger_rate's time grid as the common domain
  // This ensures all contributions are defined on the same TBP grid
  let tbp_points_sorted = integral_merger_rate.breakpoints().clone();

  // Evaluate I(t) = ∫₀ᵗ κ(t') dt' at grid points
  // I(t) represents cumulative merger rate from present (t=0) to time t
  let integral_at_node = integral_merger_rate.values().clone();

  // Compute coalescent contribution for internal node
  // An internal node with k children represents a merger event.
  // The coalescent probability density at time t is:
  //
  //   P(merger at t | k children) ∝ λ(t)^(k-1) · exp(-I(t))
  //
  // where:
  //   - λ(t) = k(t)·(k(t)-1)/(2·Tc(t)) is total merger rate
  //   - k(t) is number of concurrent lineages at time t
  //   - Tc(t) is effective population size
  //   - multiplicity = k - 1 (number of mergers = children - 1)
  //
  // In neg-log space, the formula becomes:
  //   neg_log: multiplicity · (I(t) - log(λ(t)))
  //   probability: exp(-multiplicity · (I(t) - log(λ(t))))
  //              = exp(-multiplicity · I(t)) · λ(t)^multiplicity
  //
  // Node contribution is MULTIPLIED with the product of child messages
  // during message passing.

  let n_children = node.child_edges.len() as f64;
  let multiplicity = n_children - 1.0;

  // Compute λ(t) = k(t)·(k(t)-1)/(2·Tc(t)) at grid points
  let k_vals = lineage_counts.eval_many(&tbp_points_sorted);
  let tc_vals = tc_dist.eval_many(&tbp_points_sorted)?;
  let (_, total_rate) = compute_merger_rates(&k_vals, &tc_vals);

  // Compute contribution in probability space:
  // exp(-multiplicity · (I(t) - log(λ(t))))
  // = λ(t)^multiplicity · exp(-multiplicity · I(t))
  let mut contrib_values = Array1::zeros(tbp_points_sorted.len());
  for i in 0..tbp_points_sorted.len() {
    let i_t = integral_at_node[i];
    let lambda_t = total_rate[i];
    let log_lambda_t = lambda_t.ln();
    // neg-log contribution: multiplicity · (I(t) - log(λ(t)))
    let neg_log_contrib = multiplicity * (i_t - log_lambda_t);
    // Convert to probability
    contrib_values[i] = (-neg_log_contrib).exp();
  }

  // Return contribution in TBP coordinates (as expected by tests and for consistency with Python v0)
  // TBP grid is already sorted ascending, and contrib_values are computed on that grid
  Ok(Distribution::Function(DistributionFunction::from_arrays_nonuniform(
    &tbp_points_sorted,
    &contrib_values,
  )?))
}
