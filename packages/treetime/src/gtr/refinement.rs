use crate::ancestral::marginal::{marginal_backward, update_marginal};
use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, infer_gtr_impl};
use crate::partition::traits::{HasGtr, PartitionMarginalPasses, TransitionCounting};
use eyre::Report;
use log::{debug, info, warn};
use ndarray::Array1;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};

pub fn refine_gtr_iterative<N, E, P>(
  graph: &Graph<N, E, ()>,
  partition: &Arc<RwLock<P>>,
  iterations: usize,
  fixed_pi: Option<&Array1<f64>>,
  pc: f64,
  sampling_bias_correction: Option<f64>,
  optimize_rate: bool,
) -> Result<f64, Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: TransitionCounting<N, E> + PartitionMarginalPasses<N, E> + HasGtr,
{
  let n_states = partition.read_arc().gtr().pi.len();
  let options = InferGtrOptions {
    fixed_pi: fixed_pi.cloned(),
    pc,
    ..InferGtrOptions::default()
  };

  let counts = partition.read_arc().count_transitions(graph)?;
  let result = infer_gtr_impl(&counts, &options)?;
  *partition.write_arc().gtr_mut() = build_gtr_from_inference(n_states, &result)?;
  debug!(
    "GTR refinement: initial inference, mu = {:.6}",
    partition.read_arc().gtr().mu
  );

  if optimize_rate {
    optimize_gtr_rate(graph, partition)?;
    debug!(
      "GTR refinement: initial rate optimization, mu = {:.6}",
      partition.read_arc().gtr().mu
    );
  }

  for i in 0..iterations {
    let counts = partition.read_arc().count_transitions(graph)?;
    let result = infer_gtr_impl(&counts, &options)?;
    *partition.write_arc().gtr_mut() = build_gtr_from_inference(n_states, &result)?;

    if optimize_rate {
      optimize_gtr_rate(graph, partition)?;
    }
    debug!(
      "GTR refinement: iteration {i}, mu = {:.6}",
      partition.read_arc().gtr().mu
    );
  }

  if let Some(correction) = sampling_bias_correction {
    partition.write_arc().gtr_mut().mu *= correction;
    info!(
      "Applied sampling bias correction {correction:.4}, mu = {:.6}",
      partition.read_arc().gtr().mu
    );
  }

  let partitions = std::slice::from_ref(partition);
  let log_lh = update_marginal(graph, partitions)?;

  let guard = partition.read_arc();
  let gtr = guard.gtr();
  info!(
    "GTR refinement: final log likelihood = {log_lh:.4}, mu = {:.6}, pi = {:?}",
    gtr.mu, gtr.pi
  );

  Ok(log_lh)
}

fn build_gtr_from_inference(n_states: usize, result: &InferGtrResult) -> Result<GTR, Report> {
  GTR::new(GTRParams {
    n_states,
    mu: result.mu,
    W: Some(result.W.clone()),
    pi: result.pi.clone(),
  })
}

fn optimize_gtr_rate<N, E, P>(graph: &Graph<N, E, ()>, partition: &Arc<RwLock<P>>) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + HasGtr,
{
  let old_mu = partition.read_arc().gtr().mu;
  let sqrt_old_mu = old_mu.sqrt();
  let root_key = graph.get_exactly_one_root()?.read_arc().key();

  let lo = 0.01 * sqrt_old_mu;
  let hi = 100.0 * sqrt_old_mu;

  let partitions = std::slice::from_ref(partition);

  let run_backward = |sqrt_mu: f64| -> f64 {
    {
      let mut guard = partition.write_arc();
      guard.gtr_mut().mu = sqrt_mu * sqrt_mu;
      guard.reset_node_log_likelihoods();
    }
    match marginal_backward(graph, partitions) {
      Ok(()) => -partition.read_arc().get_log_lh(root_key),
      Err(e) => {
        warn!(
          "GTR rate optimization: backward pass failed at mu={:.6}: {e}",
          sqrt_mu * sqrt_mu
        );
        f64::INFINITY
      },
    }
  };

  let cost_lo = run_backward(lo);
  let cost_mid = run_backward(sqrt_old_mu);
  let cost_hi = run_backward(hi);

  if cost_mid < cost_lo && cost_mid < cost_hi {
    let optimal_sqrt_mu = {
      let mut cost_fn = |sqrt_mu: f64| -> f64 { run_backward(sqrt_mu) };
      brent_minimize_bracketed(&mut cost_fn, lo, sqrt_old_mu, hi, cost_mid, 1e-8, 100)
    };
    run_backward(optimal_sqrt_mu);
    debug!(
      "GTR rate optimization: optimized mu = {:.6} (from {:.6})",
      partition.read_arc().gtr().mu,
      old_mu
    );
  } else {
    partition.write_arc().gtr_mut().mu = old_mu;
    debug!("GTR rate optimization: skipped (no bracket), keeping mu = {old_mu:.6}");
  }

  Ok(())
}

#[allow(clippy::many_single_char_names, clippy::float_cmp)]
fn brent_minimize_bracketed<F>(f: &mut F, xa: f64, xb: f64, xc: f64, fb: f64, tol: f64, max_iter: usize) -> f64
where
  F: FnMut(f64) -> f64,
{
  const GOLDEN: f64 = 0.381_966_011_250_105;
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
  fn test_gtr_refinement_brent_minimize_bracketed_quadratic() {
    let mut f = |x: f64| (x - 3.0) * (x - 3.0);
    let (xa, xb, xc) = (0.0, 3.5, 10.0);
    let fb = f(xb);
    let result = brent_minimize_bracketed(&mut f, xa, xb, xc, fb, 1e-10, 100);
    assert_abs_diff_eq!(result, 3.0, epsilon = 1e-8);
  }

  #[test]
  fn test_gtr_refinement_brent_minimize_bracketed_shifted_quadratic() {
    let mut f = |x: f64| (x - 0.7) * (x - 0.7) + 1.0;
    let (xa, xb, xc) = (0.0, 1.0, 5.0);
    let fb = f(xb);
    let result = brent_minimize_bracketed(&mut f, xa, xb, xc, fb, 1e-10, 100);
    assert_abs_diff_eq!(result, 0.7, epsilon = 1e-8);
  }

  #[test]
  fn test_gtr_refinement_brent_minimize_bracketed_cosine() {
    let mut f = |x: f64| x.cos();
    let (xa, xb, xc) = (2.0, 3.0, 5.0);
    let fb = f(xb);
    let result = brent_minimize_bracketed(&mut f, xa, xb, xc, fb, 1e-10, 100);
    assert_abs_diff_eq!(result, std::f64::consts::PI, epsilon = 1e-8);
  }

  #[test]
  fn test_gtr_refinement_brent_minimize_bracketed_starts_near_interior() {
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
