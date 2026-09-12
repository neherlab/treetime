#[cfg(test)]
mod __tests__;

use crate::optimize::branch_length::validate_branch_length_value;
use crate::optimize::likelihood::OptimizationMetrics;
use crate::partition::traits::PartitionOptimizeOps;
use eyre::Report;
use rayon::prelude::*;
use statrs::function::factorial::ln_factorial;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_primitives::LogLh;
use treetime_utils::{make_error, make_report};

/// Poisson indel log-likelihood contribution for one edge.
///
/// Given $k$ observed indel events on a branch of length $t$ with indel rate $\mu$
/// (indels per unit branch length), the Poisson log-likelihood and its derivatives are:
///
/// $$\ell(t) = k \ln(\mu t) - \mu t - \ln(k!)$$
/// $$\frac{d\ell}{dt} = \frac{k}{t} - \mu$$
/// $$\frac{d^2\ell}{dt^2} = -\frac{k}{t^2}$$
///
/// When $k = 0$: $\ell(t) = -\mu t$, $d\ell/dt = -\mu$, $d^2\ell/dt^2 = 0$.
/// When $k > 0$ and $t \to 0^+$: $d\ell/dt \to +\infty$, forcing the optimum away from zero.
pub fn poisson_indel_log_lh(k: usize, mu: f64, t: f64) -> Result<OptimizationMetrics, Report> {
  validate_branch_length_value(t)?;

  if mu == 0.0 {
    return Ok(if k == 0 {
      OptimizationMetrics::default()
    } else {
      OptimizationMetrics::new(LogLh::IMPOSSIBLE, 0.0, 0.0)
    });
  }

  if k == 0 {
    // Poisson(0 | mu*t) = exp(-mu*t)
    return Ok(OptimizationMetrics::new(LogLh::new(-mu * t), -mu, 0.0));
  }

  if t == 0.0 {
    return make_error!("Poisson indel likelihood requires a positive branch length when k > 0, got t={t}");
  }

  let k_f = k as f64;
  let lambda = mu * t;
  let log_lh = k_f * lambda.ln() - lambda - ln_factorial(k as u64);
  let derivative = k_f / t - mu;
  let second_derivative = -k_f / (t * t);

  Ok(OptimizationMetrics::new(
    LogLh::new(log_lh),
    derivative,
    second_derivative,
  ))
}

/// Estimate the global indel rate from the tree.
///
/// $\hat{\mu} = \frac{\sum_e k_e}{\sum_e t_e}$
///
/// where $k_e$ is the indel count on edge $e$ and $t_e$ is the branch length.
/// Returns 0 when there are no indels or total branch length is zero.
///
/// The branch length of
/// each edge is supplied by the `branch_lengths` value map (ancestral, timetree, ...).
pub fn estimate_indel_rate(
  graph: &Graph<()>,
  partitions: &[&dyn PartitionOptimizeOps],
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> f64 {
  let per_edge = graph
    .get_edges()
    .par_iter()
    .map(|edge_ref| {
      let edge_key = edge_ref.read_arc().key();
      let branch_length = branch_lengths[&edge_key].unwrap_or(0.0);
      let edge_indels = partitions.iter().map(|p| p.edge_indel_count(edge_key)).sum::<usize>();
      (edge_indels, branch_length)
    })
    .collect::<Vec<_>>();
  let total_indels = per_edge.iter().map(|(indels, _)| indels).sum::<usize>();
  let total_branch_length = per_edge.iter().map(|(_, branch_length)| branch_length).sum::<f64>();

  if total_branch_length > 0.0 && total_indels > 0 {
    total_indels as f64 / total_branch_length
  } else {
    0.0
  }
}

/// Sum the Poisson indel log-likelihood over all edges in the graph.
///
/// Uses the same `edge_indel_count()` aggregation and `poisson_indel_log_lh()`
/// evaluator as the per-edge branch-length optimizer, but evaluated at the
/// tree's current branch lengths. For an indel-bearing edge at zero branch
/// length, the Poisson log-likelihood is $-\infty$.
pub fn total_indel_log_lh(
  graph: &Graph<()>,
  partitions: &[&dyn PartitionOptimizeOps],
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  indel_rate: f64,
) -> Result<LogLh, Report> {
  graph
    .get_edges()
    .par_iter()
    .map(|edge_ref| -> Result<LogLh, Report> {
      let edge_key = edge_ref.read_arc().key();
      let branch_length = branch_lengths[&edge_key].ok_or_else(|| {
        make_report!("Cannot evaluate indel likelihood for edge {edge_key} with a missing branch length")
      })?;
      validate_branch_length_value(branch_length)?;
      let indel_count: usize = partitions.iter().map(|p| p.edge_indel_count(edge_key)).sum();
      if indel_count > 0 && branch_length <= 0.0 {
        Ok(LogLh::IMPOSSIBLE)
      } else {
        Ok(poisson_indel_log_lh(indel_count, indel_rate, branch_length)?.log_lh)
      }
    })
    .collect::<Vec<_>>()
    .into_iter()
    .sum()
}
