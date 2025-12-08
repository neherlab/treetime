use crate::commands::timetree::coalescent::integration::compute_merger_rates;
use crate::commands::timetree::coalescent::piecewise_constant_fn::PiecewiseConstantFn;
use crate::commands::timetree::coalescent::piecewise_linear_fn::PiecewiseLinearFn;
use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_formula::DistributionFormula;
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
///   Represents survival probability. Python v0 stores -I(t) in neg-log space,
///   which converts to exp(I(t)) in probability space.
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
    // Skip nodes without time distributions
    if node.payload.time_distribution.is_none() {
      return;
    }

    let contrib = if node.is_leaf {
      compute_leaf_contribution_single(&node, integral_merger_rate)
    } else {
      compute_node_contribution_single(&node, integral_merger_rate, tc_dist, lineage_counts, present_time)
    }
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

fn compute_leaf_contribution_single(
  _node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
  integral_merger_rate: &PiecewiseLinearFn,
) -> Result<Distribution, Report> {
  // Leaf nodes represent sampled lineages. The coalescent contribution encodes
  // the survival probability: the probability that this lineage existed without
  // coalescing from present (t=0) to its sampling time.
  //
  // In Python v0 with neg-log representation: -I(t) is stored, which converts to exp(I(t))
  // In Rust probability space: we need exp(-I(t)) but Python stores it as exp(I(t))
  // due to sign convention differences. Following Python exactly: exp(I(t))

  let t_min = integral_merger_rate.breakpoints()[0];
  let t_max = integral_merger_rate.breakpoints()[integral_merger_rate.breakpoints().len() - 1];

  let integral_merger_rate = Arc::new(integral_merger_rate.clone());

  let eval_fn = move |t: f64| -> eyre::Result<f64> {
    let i_t = integral_merger_rate.eval(t);
    // Python stores -I(t) in neg-log space, converts to exp(I(t)) in probability
    // Clamp to avoid overflow: exp(x) overflows for x > 700
    let i_t_clamped = i_t.min(700.0);
    // Match Python exactly: exp(I(t))
    Ok(i_t_clamped.exp())
  };

  Ok(Distribution::Formula(DistributionFormula::new(eval_fn, t_min, t_max)))
}

fn compute_node_contribution_single(
  node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
  integral_merger_rate: &PiecewiseLinearFn,
  tc_dist: &Distribution,
  lineage_counts: &PiecewiseConstantFn,
  _present_time: f64,
) -> Result<Distribution, Report> {
  // Compute coalescent contribution for internal node using exact formula evaluation.
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

  let t_min = integral_merger_rate.breakpoints()[0];
  let t_max = integral_merger_rate.breakpoints()[integral_merger_rate.breakpoints().len() - 1];

  // Clone Arc-wrapped data for use in closure
  let integral_merger_rate = Arc::new(integral_merger_rate.clone());
  let lineage_counts = Arc::new(lineage_counts.clone());
  let tc_dist = Arc::new(tc_dist.clone());

  // Create closure that evaluates coalescent contribution at any time point
  let eval_fn = move |t: f64| -> eyre::Result<f64> {
    let i_t = integral_merger_rate.eval(t);
    let k_t = lineage_counts.eval(t);
    let tc_t = tc_dist.eval(t)?;

    let (_, lambda_t) = compute_merger_rates(&Array1::from_vec(vec![k_t]), &Array1::from_vec(vec![tc_t]));
    let log_lambda_t = lambda_t[0].ln();

    // neg-log contribution: multiplicity · (I(t) - log(λ(t)))
    let neg_log_contrib = multiplicity * (i_t - log_lambda_t);
    // Convert to probability
    Ok((-neg_log_contrib).exp())
  };

  Ok(Distribution::Formula(DistributionFormula::new(eval_fn, t_min, t_max)))
}
