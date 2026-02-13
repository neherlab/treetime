use crate::commands::timetree::coalescent::integration::compute_merger_rates;
use crate::commands::timetree::coalescent::piecewise_constant_fn::PiecewiseConstantFn;
use crate::commands::timetree::coalescent::piecewise_linear_fn::PiecewiseLinearFn;
use crate::commands::timetree::timetree_traits::TimetreeNode;
use crate::distribution::distribution::{Distribution, DistributionNegLog};
use crate::distribution::distribution_formula::DistributionFormula;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use eyre::{Context, Report};
use indexmap::IndexMap;
use ndarray::Array1;
use std::sync::Arc;

/// Computes coalescent prior contributions for all nodes.
///
/// Returns map from node key to distribution that should be multiplied
/// with node's time distribution during backward pass.
///
/// Kingman coalescent contribution formula (NegLog space):
/// - For leaves: -I(t) where I(t) = ∫₀ᵗ κ(t')dt' is cumulative merger rate.
///   Probability P = exp(I(t)), so NegLog = -ln(P) = -I(t).
/// - For internal nodes with k children: multiplicity · (I(t) - log(λ(t)))
///   where λ(t) = k(t)·(k(t)-1)/(2·Tc(t)) is total merger rate
///   and m = k - 1 is multiplicity.
pub fn compute_node_contributions<N, E, D>(
  graph: &Graph<N, E, D>,
  integral_merger_rate: &PiecewiseLinearFn,
  tc_dist: &Distribution,
  lineage_counts: &PiecewiseConstantFn,
  present_time: f64,
) -> Result<IndexMap<GraphNodeKey, Arc<DistributionNegLog>>, Report>
where
  N: GraphNode + TimetreeNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  let mut contributions = IndexMap::new();

  graph.iter_breadth_first_forward(|node| {
    // Skip nodes without time distributions
    if node.payload.time_distribution().is_none() {
      return;
    }

    let contrib = if node.is_leaf {
      compute_leaf_contribution_single(integral_merger_rate)
    } else {
      compute_internal_contribution_single(
        node.child_edges.len(),
        integral_merger_rate,
        tc_dist,
        lineage_counts,
        present_time,
      )
    }
    .wrap_err_with(|| {
      let name = node.payload.name();
      let name = name.as_ref().map_or("", |n| n.as_ref());
      format!(
        "When computing coalescent contributions for node (\"{name}\") (#{})",
        node.key,
      )
    })
    .unwrap();

    contributions.insert(node.key, Arc::new(contrib));
  });

  Ok(contributions)
}

fn compute_leaf_contribution_single(integral_merger_rate: &PiecewiseLinearFn) -> Result<DistributionNegLog, Report> {
  // Leaf nodes represent sampled lineages. The coalescent contribution encodes
  // the survival probability P = exp(I(t)).
  //
  // In NegLog space: -ln(P) = -I(t).

  let t_min = integral_merger_rate.breakpoints()[0];
  let t_max = integral_merger_rate.breakpoints()[integral_merger_rate.breakpoints().len() - 1];

  let integral_merger_rate = Arc::new(integral_merger_rate.clone());

  let eval_fn = move |t: f64| -> eyre::Result<f64> {
    let i_t = integral_merger_rate.eval(t);
    // Return -I(t) for NegLog space
    Ok(-i_t)
  };

  Ok(Distribution::Formula(DistributionFormula::new(eval_fn, t_min, t_max)))
}

fn compute_internal_contribution_single(
  n_children: usize,
  integral_merger_rate: &PiecewiseLinearFn,
  tc_dist: &Distribution,
  lineage_counts: &PiecewiseConstantFn,
  _present_time: f64,
) -> Result<DistributionNegLog, Report> {
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

  let n_children = n_children as f64;
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
    Ok(neg_log_contrib)
  };

  Ok(Distribution::Formula(DistributionFormula::new(eval_fn, t_min, t_max)))
}
