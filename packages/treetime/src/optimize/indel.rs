#[cfg(test)]
mod __tests__;

use crate::optimize::branch_length::validate_branch_length_value;
use crate::optimize::likelihood::OptimizationMetrics;
use eyre::Report;
use rayon::prelude::*;
use statrs::function::factorial::ln_factorial;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_primitives::LogLh;
use treetime_utils::{make_error, make_report};

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub fn estimate_indel_rate(
  graph: &Graph,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> f64 {
  let per_edge = graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| {
      let edge_key = edge_ref.key();
      let branch_length = branch_lengths[&edge_key].unwrap_or(0.0);
      let edge_indels = indel_counts[&edge_key];
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

pub fn total_indel_log_lh(
  graph: &Graph,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  indel_rate: f64,
) -> Result<LogLh, Report> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| -> Result<LogLh, Report> {
      let edge_key = edge_ref.key();
      let branch_length = branch_lengths[&edge_key].ok_or_else(|| {
        make_report!("Cannot evaluate indel likelihood for edge {edge_key} with a missing branch length")
      })?;
      validate_branch_length_value(branch_length)?;
      let indel_count: usize = indel_counts[&edge_key];
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

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
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
