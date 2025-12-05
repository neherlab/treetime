use crate::commands::timetree::coalescent::integration::compute_merger_rates;
use crate::commands::timetree::coalescent::lineage_dynamics::PiecewiseConstant;
use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_map::distribution_map;
use crate::distribution::distribution_resample::distribution_resample;
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
  integral_merger_rate: &Distribution,
  tc_dist: &Distribution,
  lineage_counts: &PiecewiseConstant,
  present_time: f64,
) -> Result<IndexMap<GraphNodeKey, Arc<Distribution>>, Report> {
  let mut contributions = IndexMap::new();

  graph.iter_breadth_first_forward(|node| {
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
  integral_merger_rate: &Distribution,
  tc_dist: &Distribution,
  lineage_counts: &PiecewiseConstant,
  present_time: f64,
) -> Result<Distribution, Report> {
  // Get node's time distribution domain (calendar time coordinates)
  let time_points = &node
    .payload
    .time_distribution
    .as_ref()
    .expect("Node missing time distribution")
    .t();

  // Convert calendar time to time-before-present (TBP) coordinates
  // TBP increases going backward in time: TBP = present_time - calendar_time
  let tbp_points = time_points.mapv(|t| present_time - t);

  // Sort TBP points in ascending order for distribution operations
  // (time_points is sorted ascending in calendar time → tbp_points is descending)
  let mut tbp_points_sorted_vec = tbp_points.to_vec();
  tbp_points_sorted_vec.sort_by(|a, b| a.partial_cmp(b).unwrap());
  let tbp_points_sorted = Array1::from(tbp_points_sorted_vec);

  // Evaluate I(t) = ∫₀ᵗ κ(t') dt' at node time points
  // I(t) represents cumulative merger rate from present (t=0) to time t
  let integral_at_node = distribution_resample(integral_merger_rate, &tbp_points_sorted)?;

  // Compute coalescent contribution based on node type
  let contrib = if node.is_leaf {
    // LEAF NODE CONTRIBUTION
    // In neg-log space: y = -ln(P) = -I(t)
    // Converting to probability: P = exp(-y) = exp(-(-I(t))) = exp(I(t))
    //
    // This seemingly counterintuitive positive exponent arises because:
    // - I(t) represents cumulative hazard/merger rate from present to time t
    // - neg-log probabilities are used throughout
    // - The minus sign in "-I(t)" negates the usual negative log-probability convention
    // - Result: leaves contribute exp(I(t)) > 1 in probability space
    distribution_map(&integral_at_node, |i_t| i_t.exp())?
  } else {
    // INTERNAL NODE CONTRIBUTION
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

    // Compute λ(t) = k(t)·(k(t)-1)/(2·Tc(t)) at node time points
    let k_vals = lineage_counts.eval_many(&tbp_points_sorted);
    let tc_vals = tc_dist.eval_many(&tbp_points_sorted)?;
    let (_, total_rate) = compute_merger_rates(&k_vals, &tc_vals);

    // Compute contribution in probability space:
    // exp(-multiplicity · (I(t) - log(λ(t))))
    // = λ(t)^multiplicity · exp(-multiplicity · I(t))
    let mut contrib_values = Array1::zeros(tbp_points_sorted.len());
    for i in 0..tbp_points_sorted.len() {
      let i_t = integral_at_node.y()[i];
      let lambda_t = total_rate[i];
      let log_lambda_t = lambda_t.ln();
      // neg-log contribution: multiplicity · (I(t) - log(λ(t)))
      let neg_log_contrib = multiplicity * (i_t - log_lambda_t);
      // Convert to probability
      contrib_values[i] = (-neg_log_contrib).exp();
    }

    Distribution::function(tbp_points_sorted, contrib_values)?
  };

  // Resample contribution back to original unsorted tbp_points
  // (which correspond to time_points in calendar coordinates)
  let contrib_values = contrib.eval_many(&tbp_points)?;
  Distribution::function(time_points.clone(), contrib_values)
}
