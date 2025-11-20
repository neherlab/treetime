use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_map::distribution_map;
use crate::distribution::distribution_multiplication::distribution_multiplication;
use crate::distribution::distribution_resample::distribution_resample;
use crate::graph::graph::GraphNodeForward;
use crate::graph::node::GraphNodeKey;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use eyre::{Context, Report};
use indexmap::IndexMap;
use ndarray::{Array1, s};
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_utils::make_error;

/// Computes Kingman coalescent prior contributions for all nodes in the phylogenetic tree.
///
/// Returns distributions that encode coalescent likelihood contributions for each node,
/// to be multiplied with node time distributions during backward pass optimization.
///
/// # Meaning
///
/// The Kingman coalescent model provides a probabilistic framework for assessing whether
/// the timing and structure of a phylogenetic tree are consistent with a given population
/// history. Not all tree topologies and divergence times are equally likely - trees with
/// many lineages coalescing simultaneously are less probable than gradual coalescence,
/// especially in large populations.
///
/// This function computes how likely each node's divergence time is under the coalescent
/// model, given the population size history Tc(t) and the number of concurrent lineages k(t).
/// These likelihood contributions act as priors that guide the time tree inference toward
/// more biologically plausible configurations by penalizing unlikely coalescence patterns.
///
/// # Notation and Terms
///
/// - `t` - time (negative values for past, zero at present)
/// - `k(t)` - number of concurrent lineages at time t
/// - `Tc(t)` - coalescence time scale (effective population size) at time t, controls merger rate
/// - `κ(t)` - branch merger rate: rate for one branch to merge with any other = (k(t)-1)/(2*Tc(t))
/// - `λ(t)` - total merger rate: rate for any merger to occur = k(t)*(k(t)-1)/(2*Tc(t))
/// - `I(t)` - cumulative merger rate: I(t) = ∫₀ᵗ κ(t') dt'
///
/// # Kingman Coalescent Probability Density
///
/// For a node at time t with m children (m=1 for leaves, m≥2 for internal nodes):
///
/// - Leaf nodes (m=1): P(t) ∝ exp(-I(t))
///   - Probability of no merger from present to time t
///
/// - Internal nodes (m≥2): P(t) ∝ λ(t)^(m-1) · exp(-I(t))
///   - λ(t)^(m-1): probability density of m-way merger at time t
///   - exp(-I(t)): probability of no merger before time t
///
/// # Returns
///
/// Map from node keys to distributions representing coalescent prior contributions.
/// Each distribution should be multiplied with the node's time distribution.
pub fn compute_coalescent_contributions(
  graph: &GraphAncestral,
  tc: &Distribution,
) -> Result<IndexMap<GraphNodeKey, Arc<Distribution>>, Report> {
  let (present_time, events_calendar) = collect_tree_events(graph)?;
  let events_tbp: Vec<_> = events_calendar
    .iter()
    .map(|(t, delta)| (present_time - *t, *delta))
    .collect();

  let lineage_counts = compute_lineage_count_distribution(&events_tbp)?;
  let integral_merger_rate = compute_integral_merger_rate(tc, &lineage_counts)?;
  compute_node_contributions(graph, &integral_merger_rate, tc, &lineage_counts, present_time)
}

/// Collects tree merger events as (time, delta_branches) tuples.
///
/// Returns events sorted by decreasing time (past to present).
/// delta_branches: +1 for leaf nodes, -(k-1) for internal nodes with k children.
fn collect_tree_events(graph: &GraphAncestral) -> Result<(f64, Vec<(f64, i32)>), Report> {
  let mut max_time = f64::NEG_INFINITY;
  let mut events = Vec::new();

  graph.iter_breadth_first_forward(|node| {
    if let Some(time_dist) = &node.payload.time_distribution {
      if let Some(t) = time_dist.likely_time() {
        if t > max_time {
          max_time = t;
        }
        let num_children = node.child_edges.len();

        if num_children == 0 {
          events.push((t, 1));
        } else {
          events.push((t, -((num_children as i32) - 1)));
        }
      }
    }
  });

  if events.is_empty() {
    return make_error!("No tree events found");
  }

  if !max_time.is_finite() {
    return make_error!("Cannot determine present time for coalescent events");
  }

  events.sort_by_key(|x| OrderedFloat(x.0));

  Ok((max_time, events))
}

/// Computes k(t) distribution from tree events.
///
/// k(t) is the number of concurrent lineages at time t.
/// The function is piecewise constant, stepping at each merger event.
/// Events must be sorted by decreasing time (past to present).
fn compute_lineage_count_distribution(events: &[(f64, i32)]) -> Result<Distribution, Report> {
  if events.is_empty() {
    return make_error!("Cannot build lineage count interpolator from empty events");
  }

  let mut aggregated_events = BTreeMap::new();
  for &(time, delta) in events {
    *aggregated_events.entry(OrderedFloat(time)).or_insert(0) += delta;
  }

  let mut cumulative_events: Vec<(f64, i32)> = Vec::with_capacity(aggregated_events.len());
  let mut current_count = 0_i32;
  for (time, delta) in aggregated_events.into_iter().rev() {
    current_count += delta;
    cumulative_events.push((time.into_inner(), current_count));
  }
  cumulative_events.reverse();

  let data_t_min = cumulative_events.first().unwrap().0;
  let data_t_max = cumulative_events.last().unwrap().0;
  let data_range = data_t_max - data_t_min;
  let margin = f64::max(data_range * 0.1, 1.0);

  let t_min = data_t_min - margin;
  let t_max = data_t_max + margin;
  let n_points = 10000;
  let dx = (t_max - t_min) / (n_points - 1) as f64;

  let mut y_vals = Array1::zeros(n_points);
  let mut prev_grid_idx = 0;

  for (event_time, count) in cumulative_events {
    let next_grid_idx = (((event_time - t_min) / dx).ceil() as usize).min(n_points);

    if prev_grid_idx < next_grid_idx {
      y_vals.slice_mut(s![prev_grid_idx..next_grid_idx]).fill(count as f64);
    }

    prev_grid_idx = next_grid_idx;
  }

  let func = DistributionFunction::from_range_values((t_min, t_max), y_vals)?;
  Ok(Distribution::Function(func))
}

/// Computes I(t) = ∫₀ᵗ κ(t') dt' via trapezoidal integration.
///
/// This integral represents the expected number of merger events experienced by a branch.
/// Uses the exact time points from lineage_counts where the function has discontinuities.
fn compute_integral_merger_rate(tc_dist: &Distribution, lineage_counts: &Distribution) -> Result<Distribution, Report> {
  let tvals = lineage_counts.t();
  if tvals.len() < 2 {
    return make_error!("lineage count distribution must have at least 2 points");
  }

  let k_vals = lineage_counts.eval_many(&tvals)?;
  let tc_vals = tc_dist.eval_many(&tvals)?;
  let k_clamped = k_vals.mapv(|x| f64::max(0.5, x - 1.0));
  let branch_rates = 0.5 * &k_clamped / &tc_vals;

  let n_points = tvals.len();
  let mut integral_values = Array1::zeros(n_points);
  for i in 1..n_points {
    let dt = tvals[i] - tvals[i - 1];
    let avg_rate = 0.5 * (branch_rates[i - 1] + branch_rates[i]);
    integral_values[i] = integral_values[i - 1] + dt * avg_rate;
  }

  Distribution::function(tvals.to_owned(), integral_values)
}

/// Computes coalescent prior contributions for all nodes.
///
/// Returns map from node key to distribution that should be multiplied
/// with node's time distribution during backward pass.
///
/// Kingman coalescent probability density:
/// - For leaves: P(t) ∝ exp(-I(t)) where I(t) = ∫κ(t')dt' is cumulative merger rate
/// - For internal nodes with k children: P(t) ∝ λ(t)^(k-1) · exp(-I(t))
///   where λ(t) = k(k-1)/(2Tc) is total merger rate
fn compute_node_contributions(
  graph: &GraphAncestral,
  integral_merger_rate: &Distribution,
  tc_dist: &Distribution,
  lineage_counts: &Distribution,
  present_time: f64,
) -> Result<IndexMap<GraphNodeKey, Arc<Distribution>>, Report> {
  let mut contributions = IndexMap::new();

  graph.iter_breadth_first_forward(|node| {
    let contrib = compute_node_contributions_single(&node, integral_merger_rate, tc_dist, lineage_counts, present_time)
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

fn compute_node_contributions_single(
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

  let contrib = distribution_resample(&contrib, time_points)?;

  Ok(contrib)
}

/// Computes branch merger rate κ(t) and total merger rate λ(t).
///
/// κ(t) = (k(t)-1)/(2*Tc(t)) - rate for one branch to merge with any other
/// λ(t) = k(t)*(k(t)-1)/(2*Tc(t)) - rate for any branch to merge with any other
///
/// k is clamped to ensure positive rates even in edge cases.
fn compute_merger_rates(k: &Array1<f64>, tc: &Array1<f64>) -> (Array1<f64>, Array1<f64>) {
  let k_clamped = k.mapv(|x| f64::max(0.5, x - 1.0));
  let branch_rate = 0.5 * &k_clamped / tc;
  let total_rate = 0.5 * &k_clamped * (&k_clamped + 1.0) / tc;
  (branch_rate, total_rate)
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;
  use ndarray::array;

  #[test]
  fn test_lineage_count_simple_tree() -> Result<(), Report> {
    let events = vec![(10.0, 1), (10.0, 1), (5.0, -1), (0.0, 1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    // Check the actual grid values directly (no interpolation)
    let t_grid = lineage_counts.t();
    let y_grid = lineage_counts.y();

    // The grid should sample the step function densely
    // We expect roughly constant values in each region
    // Just verify it's reasonable by checking a few properties
    assert!(t_grid.len() == 10000, "Expected 10000 grid points");

    // Find grid points in different regions and check their values
    let idx_future = t_grid.iter().position(|&t| t > 10.5).unwrap();
    let idx_high = t_grid.iter().position(|&t| t > 7.0 && t < 9.0).unwrap();
    let idx_mid = t_grid.iter().position(|&t| t > 2.0 && t < 4.0).unwrap();
    let idx_past = t_grid.iter().position(|&t| t < -0.5).unwrap();

    assert_abs_diff_eq!(y_grid[idx_future], 0.0, epsilon = 0.1);
    assert_abs_diff_eq!(y_grid[idx_high], 2.0, epsilon = 0.1);
    assert_abs_diff_eq!(y_grid[idx_mid], 1.0, epsilon = 0.1);
    assert_abs_diff_eq!(y_grid[idx_past], 2.0, epsilon = 0.1);

    Ok(())
  }

  #[test]
  fn test_merger_rates() -> Result<(), Report> {
    let k = array![2.0, 3.0, 4.0];
    let tc = array![0.001, 0.002, 0.003];

    let (branch_rate, total_rate) = compute_merger_rates(&k, &tc);

    assert_abs_diff_eq!(branch_rate[0], 0.5 * 1.0 / 0.001, epsilon = 1e-10);
    assert_abs_diff_eq!(branch_rate[1], 0.5 * 2.0 / 0.002, epsilon = 1e-10);
    assert_abs_diff_eq!(branch_rate[2], 0.5 * 3.0 / 0.003, epsilon = 1e-10);

    assert_abs_diff_eq!(total_rate[0], 0.5 * 1.0 * 2.0 / 0.001, epsilon = 1e-10);
    assert_abs_diff_eq!(total_rate[1], 0.5 * 2.0 * 3.0 / 0.002, epsilon = 1e-10);
    assert_abs_diff_eq!(total_rate[2], 0.5 * 3.0 * 4.0 / 0.003, epsilon = 1e-10);

    Ok(())
  }
}
