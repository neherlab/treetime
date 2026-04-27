use crate::commands::mugration::discrete_marginal::{discrete_marginal_backward, run_discrete_marginal};
use crate::constants::SUPERTINY_NUMBER;
use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::infer_gtr::common::{
  InferGtrOptions, InferGtrResult, MutationCounts, infer_gtr_impl, is_profile_informative,
};
use crate::gtr::infer_gtr::dense::{accumulate_mutation_counts, get_branch_mutation_matrix};
use crate::make_report;
use crate::representation::partition::discrete::PartitionDiscrete;
use crate::representation::partition::traits::HasLogLh;
use eyre::Report;
use log::{debug, info, warn};
use ndarray::{Array1, Array2};
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_utils::array::ndarray::argmax_first;

/// Iteratively refine the GTR model for discrete trait reconstruction.
///
/// After an initial forward-backward pass, alternates between:
/// 1. Counting transitions from marginal profiles → re-estimating W (and pi if not fixed)
/// 2. Optimizing the scalar rate mu via Brent minimization
///
/// Follows v0's `reconstruct_discrete_traits()` workflow: initial reconstruction with GTR
/// inference, then `iterations` rounds of `infer_gtr()` + `optimize_gtr_rate()`, optionally
/// applying sampling bias correction, and a final reconstruction with the refined model.
pub fn refine_gtr_iterative<N, E>(
  graph: &Graph<N, E, ()>,
  partition: &mut PartitionDiscrete,
  iterations: usize,
  fixed_pi: Option<&Array1<f64>>,
  pc: f64,
  sampling_bias_correction: Option<f64>,
) -> Result<f64, Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let n_states = partition.n_states();
  let options = InferGtrOptions {
    fixed_pi: fixed_pi.cloned(),
    pc,
    ..InferGtrOptions::default()
  };

  // Initial GTR inference from the first reconstruction
  let counts = count_transitions_discrete(graph, partition)?;
  let result = infer_gtr_impl(&counts, &options)?;
  partition.gtr = build_gtr_from_inference(n_states, &result)?;
  debug!(
    "Mugration GTR refinement: initial inference, mu = {:.6}",
    partition.gtr.mu
  );

  // Optimize overall rate
  optimize_gtr_rate(graph, partition)?;
  debug!(
    "Mugration GTR refinement: initial rate optimization, mu = {:.6}",
    partition.gtr.mu
  );

  // Iterative refinement: re-estimate W from current profiles, then optimize mu
  for i in 0..iterations {
    let counts = count_transitions_discrete(graph, partition)?;
    let result = infer_gtr_impl(&counts, &options)?;
    partition.gtr = build_gtr_from_inference(n_states, &result)?;

    optimize_gtr_rate(graph, partition)?;
    debug!("Mugration GTR refinement: iteration {i}, mu = {:.6}", partition.gtr.mu);
  }

  // Apply sampling bias correction before final reconstruction
  if let Some(correction) = sampling_bias_correction {
    partition.gtr.mu *= correction;
    info!(
      "Mugration: applied sampling bias correction {correction:.4}, mu = {:.6}",
      partition.gtr.mu
    );
  }

  // Final reconstruction with refined model
  let log_lh = run_discrete_marginal(graph, partition)?;
  info!(
    "Mugration GTR refinement: final log likelihood = {log_lh:.4}, mu = {:.6}, pi = {:?}",
    partition.gtr.mu, partition.gtr.pi
  );

  Ok(log_lh)
}

/// Count state transitions from marginal profiles in a discrete partition.
///
/// For each edge, computes the joint parent-child probability matrix
/// `M[i,j] = P(child=i, parent=j | data)` from the marginal profiles
/// (subtree likelihood at child × outgroup likelihood at parent × transition probability),
/// then accumulates expected transition counts and dwell times.
fn count_transitions_discrete<N, E>(
  graph: &Graph<N, E, ()>,
  partition: &PartitionDiscrete,
) -> Result<MutationCounts, Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let n_states = partition.n_states();
  let mut nij = Array2::zeros((n_states, n_states));
  let mut Ti = Array1::zeros(n_states);

  for edge in graph.get_edges() {
    let edge_arc = edge.read_arc();
    let branch_length = partition.effective_branch_length(edge_arc.payload().read_arc().branch_length().unwrap_or(0.0));
    let edge_key = edge_arc.key();

    let edge_data = &partition.edges[&edge_key];

    // Reshape 1D profiles (n_states,) to 2D (1, n_states) for get_branch_mutation_matrix
    let msg_to_child = edge_data
      .msg_to_child
      .clone()
      .into_shape_with_order((1, n_states))
      .map_err(|e| make_report!("reshape msg_to_child: {e}"))?;
    let msg_to_parent = edge_data
      .msg_to_parent
      .clone()
      .into_shape_with_order((1, n_states))
      .map_err(|e| make_report!("reshape msg_to_parent: {e}"))?;

    let exp_qt = partition.gtr.expQt(branch_length) + SUPERTINY_NUMBER;
    let mut_stack = get_branch_mutation_matrix(&msg_to_child, &msg_to_parent, &exp_qt);
    accumulate_mutation_counts(&mut_stack, branch_length, &mut nij, &mut Ti);
  }

  // Root state: one-hot from argmax of root profile, skipping near-uniform profiles.
  // Near-uniform root profiles carry no phylogenetic signal; converting them to one-hot
  // injects state-order bias into pi estimation. Matches the dense GTR inference path.
  let root = graph.get_exactly_one_root()?;
  let root_key = root.read_arc().key();
  let root_profile = &partition.nodes[&root_key].profile;
  let mut root_state = Array1::zeros(n_states);
  if is_profile_informative(root_profile.view(), n_states) {
    if let Some(root_idx) = argmax_first(&root_profile.view()) {
      root_state[root_idx] = 1.0;
    }
  }

  // Zero diagonal (no-change events excluded)
  nij.diag_mut().fill(0.0);

  Ok(MutationCounts { nij, Ti, root_state })
}

fn build_gtr_from_inference(n_states: usize, result: &InferGtrResult) -> Result<GTR, Report> {
  GTR::new(GTRParams {
    n_states,
    mu: result.mu,
    W: Some(result.W.clone()),
    pi: result.pi.clone(),
  })
}

/// Optimize the overall GTR rate (mu) via Brent's method.
///
/// Minimizes negative log-likelihood over sqrt(mu) to ensure mu > 0.
/// Checks for a valid downhill bracket before optimizing. If the cost function
/// is monotonic (no interior minimum), keeps the current mu unchanged,
/// matching v0's fallback behavior when scipy's optimization fails.
fn optimize_gtr_rate<N, E>(graph: &Graph<N, E, ()>, partition: &mut PartitionDiscrete) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let old_mu = partition.gtr.mu;
  let sqrt_old_mu = old_mu.sqrt();
  let root_key = graph.get_exactly_one_root()?.read_arc().key();

  let lo = 0.01 * sqrt_old_mu;
  let hi = 100.0 * sqrt_old_mu;

  let run_backward = |partition: &mut PartitionDiscrete, sqrt_mu: f64| -> f64 {
    partition.gtr.mu = sqrt_mu * sqrt_mu;
    for node_data in partition.nodes.values_mut() {
      node_data.log_lh = 0.0;
    }
    match discrete_marginal_backward(graph, partition) {
      Ok(()) => -partition.get_log_lh(root_key),
      Err(e) => {
        warn!("Mugration: backward pass failed at mu={:.6}: {e}", sqrt_mu * sqrt_mu);
        f64::INFINITY
      },
    }
  };

  // Check for valid downhill bracket: f(mid) < f(lo) AND f(mid) < f(hi).
  // For mugration, -log_lh is often monotonically decreasing toward high mu
  // (converging to stationary distribution), so no interior minimum exists.
  // v0's scipy raises an error in this case and keeps old mu.
  let cost_lo = run_backward(partition, lo);
  let cost_mid = run_backward(partition, sqrt_old_mu);
  let cost_hi = run_backward(partition, hi);

  if cost_mid < cost_lo && cost_mid < cost_hi {
    let optimal_sqrt_mu = {
      let mut cost_fn = |sqrt_mu: f64| -> f64 { run_backward(partition, sqrt_mu) };
      brent_minimize_bracketed(&mut cost_fn, lo, sqrt_old_mu, hi, cost_mid, 1e-8, 100)
    };
    run_backward(partition, optimal_sqrt_mu);
    debug!(
      "Mugration: optimized rate mu = {:.6} (from {:.6})",
      partition.gtr.mu, old_mu
    );
  } else {
    // No valid bracket: restore mu but leave backward profiles at the last bracket
    // evaluation (hi). v0's scipy also leaves profiles at its last evaluation before
    // failing, then restores mu without re-running backward.
    partition.gtr.mu = old_mu;
    debug!("Mugration: rate optimization skipped (no bracket), keeping mu = {old_mu:.6}");
  }

  Ok(())
}

/// Brent's method for 1D minimization with a known bracket triple `(xa, xb, xc)`.
///
/// Matches SciPy's `Brent.optimize()` initialization: starts from the interior
/// point `xb` where `f(xb) < f(xa)` and `f(xb) < f(xc)`, preserving the
/// three-point bracket semantics that plain interval-based Brent discards.
#[allow(clippy::many_single_char_names, clippy::float_cmp)]
fn brent_minimize_bracketed<F>(f: &mut F, xa: f64, xb: f64, xc: f64, fb: f64, tol: f64, max_iter: usize) -> f64
where
  F: FnMut(f64) -> f64,
{
  const GOLDEN: f64 = 0.381_966_011_250_105; // (3 - sqrt(5)) / 2
  const ZEPS: f64 = 1e-11;

  debug_assert!(
    xb > xa.min(xc) && xb < xa.max(xc),
    "xb={xb} must be interior to [{}, {}]",
    xa.min(xc),
    xa.max(xc)
  );

  let (mut a, mut b) = if xa < xc { (xa, xc) } else { (xc, xa) };
  let mut x = xb;
  let mut w = xb;
  let mut v = xb;
  let mut fx = fb;
  let mut fw = fb;
  let mut fv = fb;
  let mut e: f64 = 0.0;
  let mut d: f64 = 0.0;

  for _ in 0..max_iter {
    let m = 0.5 * (a + b);
    let tol1 = tol * x.abs() + ZEPS;
    let tol2 = 2.0 * tol1;

    if (x - m).abs() <= tol2 - 0.5 * (b - a) {
      return x;
    }

    let mut use_golden = true;
    if e.abs() > tol1 {
      let r = (x - w) * (fx - fv);
      let q = (x - v) * (fx - fw);
      let mut p = (x - v) * q - (x - w) * r;
      let mut q = 2.0 * (q - r);
      if q > 0.0 {
        p = -p;
      } else {
        q = -q;
      }
      let old_e = e;
      e = d;

      if p.abs() < (0.5 * q * old_e).abs() && p > q * (a - x) && p < q * (b - x) {
        d = p / q;
        let u = x + d;
        if (u - a) < tol2 || (b - u) < tol2 {
          d = if m > x { tol1 } else { -tol1 };
        }
        use_golden = false;
      }
    }

    if use_golden {
      e = if x < m { b - x } else { a - x };
      d = GOLDEN * e;
    }

    let u = if d.abs() >= tol1 {
      x + d
    } else {
      x + if d > 0.0 { tol1 } else { -tol1 }
    };
    let fu = f(u);

    if fu <= fx {
      if u < x {
        b = x;
      } else {
        a = x;
      }
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    } else {
      if u < x {
        a = u;
      } else {
        b = u;
      }
      if fu <= fw || w == x {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      } else if fu <= fv || v == x || v == w {
        v = u;
        fv = fu;
      }
    }
  }

  x
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;

  #[test]
  fn test_brent_minimize_bracketed_quadratic() {
    let mut f = |x: f64| (x - 3.0) * (x - 3.0);
    let (xa, xb, xc) = (0.0, 3.5, 10.0);
    let fb = f(xb);
    let result = brent_minimize_bracketed(&mut f, xa, xb, xc, fb, 1e-10, 100);
    assert_abs_diff_eq!(result, 3.0, epsilon = 1e-8);
  }

  #[test]
  fn test_brent_minimize_bracketed_shifted_quadratic() {
    let mut f = |x: f64| (x - 0.7) * (x - 0.7) + 1.0;
    let (xa, xb, xc) = (0.0, 1.0, 5.0);
    let fb = f(xb);
    let result = brent_minimize_bracketed(&mut f, xa, xb, xc, fb, 1e-10, 100);
    assert_abs_diff_eq!(result, 0.7, epsilon = 1e-8);
  }

  #[test]
  fn test_brent_minimize_bracketed_cosine() {
    let mut f = |x: f64| x.cos();
    let (xa, xb, xc) = (2.0, 3.0, 5.0);
    let fb = f(xb);
    let result = brent_minimize_bracketed(&mut f, xa, xb, xc, fb, 1e-10, 100);
    assert_abs_diff_eq!(result, std::f64::consts::PI, epsilon = 1e-8);
  }

  #[test]
  fn test_brent_minimize_bracketed_starts_near_interior() {
    let call_count = std::cell::Cell::new(0_u32);
    let mut f = |x: f64| -> f64 {
      call_count.set(call_count.get() + 1);
      (x - 5.0) * (x - 5.0)
    };
    let (xa, xb, xc) = (0.0, 4.9, 10.0);
    let fb = f(xb);
    call_count.set(0);
    let result = brent_minimize_bracketed(&mut f, xa, xb, xc, fb, 1e-10, 100);
    assert_abs_diff_eq!(result, 5.0, epsilon = 1e-8);
    let count = call_count.get();
    assert!(
      count < 20,
      "Starting near minimum should converge quickly, took {count} evaluations"
    );
  }
}
