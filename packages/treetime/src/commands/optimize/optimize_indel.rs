use crate::commands::optimize::optimize_unified::OptimizationMetrics;
use crate::commands::optimize::partition_ops::PartitionOptimizeOps;
use crate::representation::payload::ancestral::GraphAncestral;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::HasBranchLength;

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
pub fn poisson_indel_log_lh(k: usize, mu: f64, t: f64) -> OptimizationMetrics {
  if mu == 0.0 {
    return OptimizationMetrics::default();
  }

  if k == 0 {
    // Poisson(0 | mu*t) = exp(-mu*t)
    return OptimizationMetrics::new(-mu * t, -mu, 0.0);
  }

  let k_f = k as f64;
  let lambda = mu * t;
  let log_lh = k_f * lambda.ln() - lambda - ln_factorial(k);
  let derivative = k_f / t - mu;
  let second_derivative = -k_f / (t * t);

  OptimizationMetrics::new(log_lh, derivative, second_derivative)
}

/// Estimate the global indel rate from the tree.
///
/// $\hat{\mu} = \frac{\sum_e k_e}{\sum_e t_e}$
///
/// where $k_e$ is the indel count on edge $e$ and $t_e$ is the branch length.
/// Returns 0 when there are no indels or total branch length is zero.
pub fn estimate_indel_rate<P>(graph: &GraphAncestral, partitions: &[Arc<RwLock<P>>]) -> f64
where
  P: PartitionOptimizeOps + ?Sized,
{
  let mut total_indels: usize = 0;
  let mut total_branch_length: f64 = 0.0;

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.read_arc().key();
    let branch_length = edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0);

    let edge_indels: usize = partitions.iter().map(|p| p.read_arc().edge_indel_count(edge_key)).sum();

    total_indels += edge_indels;
    total_branch_length += branch_length;
  }

  if total_branch_length > 0.0 && total_indels > 0 {
    total_indels as f64 / total_branch_length
  } else {
    0.0
  }
}

/// $\ln(k!)$ for non-negative integer $k$.
fn ln_factorial(k: usize) -> f64 {
  (1..=k).fold(0.0, |acc, i| acc + (i as f64).ln())
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;

  #[test]
  fn test_optimize_indel_poisson_zero_rate() {
    let metrics = poisson_indel_log_lh(3, 0.0, 0.1);
    assert_abs_diff_eq!(metrics.log_lh, 0.0, epsilon = 1e-15);
    assert_abs_diff_eq!(metrics.derivative, 0.0, epsilon = 1e-15);
    assert_abs_diff_eq!(metrics.second_derivative, 0.0, epsilon = 1e-15);
  }

  #[test]
  fn test_optimize_indel_poisson_zero_indels() {
    let mu = 5.0;
    let t = 0.1;
    let metrics = poisson_indel_log_lh(0, mu, t);
    assert_abs_diff_eq!(metrics.log_lh, -mu * t, epsilon = 1e-15);
    assert_abs_diff_eq!(metrics.derivative, -mu, epsilon = 1e-15);
    assert_abs_diff_eq!(metrics.second_derivative, 0.0, epsilon = 1e-15);
  }

  #[test]
  fn test_optimize_indel_poisson_log_lh_value() {
    // k=2, mu=10, t=0.1 => lambda=1.0
    // log P(2|1.0) = 2*ln(1) - 1 - ln(2!) = 0 - 1 - ln(2) = -1 - 0.6931...
    let metrics = poisson_indel_log_lh(2, 10.0, 0.1);
    let expected_log_lh = -1.0 - 2.0_f64.ln();
    assert_abs_diff_eq!(metrics.log_lh, expected_log_lh, epsilon = 1e-14);
  }

  #[test]
  fn test_optimize_indel_poisson_derivative() {
    let k = 3;
    let mu = 5.0;
    let t = 0.2;
    let metrics = poisson_indel_log_lh(k, mu, t);
    // d/dt = k/t - mu = 3/0.2 - 5 = 15 - 5 = 10
    assert_abs_diff_eq!(metrics.derivative, 10.0, epsilon = 1e-13);
    // d2/dt2 = -k/t^2 = -3/0.04 = -75
    // 0.2 is not exactly representable in float64, so t*t has rounding error
    assert_abs_diff_eq!(metrics.second_derivative, -75.0, epsilon = 1e-12);
  }

  #[test]
  fn test_optimize_indel_poisson_mle_at_optimum() {
    // At the MLE t = k/mu, the derivative should be zero
    let k = 5;
    let mu = 10.0;
    let t_mle = k as f64 / mu; // 0.5
    let metrics = poisson_indel_log_lh(k, mu, t_mle);
    assert_abs_diff_eq!(metrics.derivative, 0.0, epsilon = 1e-14);
  }

  #[test]
  fn test_optimize_indel_poisson_derivative_positive_near_zero() {
    // For k > 0, derivative should be large positive near t=0
    let metrics = poisson_indel_log_lh(1, 5.0, 1e-6);
    assert!(metrics.derivative > 1e5);
  }

  #[test]
  fn test_optimize_indel_poisson_second_derivative_negative() {
    // Second derivative is always negative when k > 0 (log-concave)
    let metrics = poisson_indel_log_lh(3, 5.0, 0.5);
    assert!(metrics.second_derivative < 0.0);
  }

  #[test]
  fn test_optimize_indel_ln_factorial() {
    assert_abs_diff_eq!(ln_factorial(0), 0.0, epsilon = 1e-15);
    assert_abs_diff_eq!(ln_factorial(1), 0.0, epsilon = 1e-15);
    assert_abs_diff_eq!(ln_factorial(2), 2.0_f64.ln(), epsilon = 1e-15);
    assert_abs_diff_eq!(ln_factorial(5), 120.0_f64.ln(), epsilon = 1e-13);
    assert_abs_diff_eq!(ln_factorial(10), 3628800.0_f64.ln(), epsilon = 1e-11);
  }
}
