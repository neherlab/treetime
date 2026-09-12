use crate::ancestral::marginal::{marginal_backward, marginal_update, profile_branch_lengths};
use crate::gtr::brent_bracketed::BrentBracketed;
use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, infer_gtr_impl};
use crate::make_internal_report;
use crate::partition::traits::{HasGtr, PartitionMarginalPasses, TransitionCounting};
use argmin::core::{CostFunction, Error, Executor};
use eyre::Report;
use log::{debug, info, warn};
use ndarray::Array1;
use std::cell::RefCell;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_primitives::LogLh;

pub fn refine_gtr_iterative<P>(
  graph: &Graph,
  partition: &RefCell<P>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  iterations: usize,
  fixed_pi: Option<&Array1<f64>>,
  pc: f64,
  sampling_bias_correction: Option<f64>,
  optimize_rate: bool,
) -> Result<LogLh, Report>
where
  P: TransitionCounting + PartitionMarginalPasses + HasGtr,
{
  let n_states = partition.borrow().gtr().pi.len();
  let options = InferGtrOptions {
    fixed_pi: fixed_pi.cloned(),
    pc,
    ..InferGtrOptions::default()
  };

  let counts = partition.borrow().count_transitions(graph, branch_lengths)?;
  let result = infer_gtr_impl(&counts, &options)?;
  *partition.borrow_mut().gtr_mut() = build_gtr_from_inference(n_states, &result)?;
  debug!(
    "GTR refinement: initial inference, mu = {:.6}",
    partition.borrow().gtr().mu
  );

  if optimize_rate {
    optimize_gtr_rate(graph, partition, branch_lengths)?;
    debug!(
      "GTR refinement: initial rate optimization, mu = {:.6}",
      partition.borrow().gtr().mu
    );
  }

  for i in 0..iterations {
    let counts = partition.borrow().count_transitions(graph, branch_lengths)?;
    let result = infer_gtr_impl(&counts, &options)?;
    *partition.borrow_mut().gtr_mut() = build_gtr_from_inference(n_states, &result)?;

    if optimize_rate {
      optimize_gtr_rate(graph, partition, branch_lengths)?;
    }
    debug!("GTR refinement: iteration {i}, mu = {:.6}", partition.borrow().gtr().mu);
  }

  if let Some(correction) = sampling_bias_correction {
    partition.borrow_mut().gtr_mut().mu *= correction;
    info!(
      "Applied sampling bias correction {correction:.4}, mu = {:.6}",
      partition.borrow().gtr().mu
    );
  }

  let log_lh = marginal_update(
    graph,
    &profile_branch_lengths(branch_lengths),
    std::slice::from_mut(&mut *partition.borrow_mut()),
  )?;

  let guard = partition.borrow();
  let gtr = guard.gtr();
  info!(
    "GTR refinement: final log likelihood = {:.4}, mu = {:.6}, pi = {:?}",
    log_lh.value(),
    gtr.mu,
    gtr.pi
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

fn optimize_gtr_rate<P>(
  graph: &Graph,
  partition: &RefCell<P>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report>
where
  P: PartitionMarginalPasses + HasGtr,
{
  let old_mu = partition.borrow().gtr().mu;
  let sqrt_old_mu = old_mu.sqrt();
  let root_key = graph.get_exactly_one_root()?.read_arc().key();

  let lo = 0.01 * sqrt_old_mu;
  let hi = 100.0 * sqrt_old_mu;

  // Branch lengths stay fixed while only the substitution rate is optimized, so derive the profile
  // map once and route the same map into every backward pass the Brent search evaluates.
  let branch_lengths = profile_branch_lengths(branch_lengths);

  let cost_fn = GtrRateCostFn {
    graph,
    partition,
    branch_lengths,
    root_key,
  };

  let cost_lo = cost_fn.neg_log_lh(lo);
  let cost_mid = cost_fn.neg_log_lh(sqrt_old_mu);
  let cost_hi = cost_fn.neg_log_lh(hi);

  if cost_mid < cost_lo && cost_mid < cost_hi {
    // Seed Brent at the interior estimate `sqrt_old_mu`, matching v0's
    // scipy.optimize.brent with a length-3 bracket. argmin's BrentOpt instead
    // seeds at the golden-section point of [lo, hi], which converges to a
    // different optimum on flat likelihood surfaces. max_iters matches scipy's
    // default (500).
    let solver = BrentBracketed::new(lo, sqrt_old_mu, hi);
    let res = Executor::new(&cost_fn, solver)
      .configure(|cfg| cfg.max_iters(500))
      .run()
      .map_err(|e| make_internal_report!("GTR rate optimization: argmin BrentBracketed failed: {e}"))?;

    let optimal_sqrt_mu = res
      .state()
      .best_param
      .ok_or_else(|| make_internal_report!("GTR rate optimization: solver succeeded but reported no best_param"))?;

    cost_fn.neg_log_lh(optimal_sqrt_mu);
    debug!(
      "GTR rate optimization: optimized mu = {:.6} (from {:.6})",
      partition.borrow().gtr().mu,
      old_mu
    );
  } else {
    partition.borrow_mut().gtr_mut().mu = old_mu;
    debug!("GTR rate optimization: skipped (no bracket), keeping mu = {old_mu:.6}");
  }

  Ok(())
}

struct GtrRateCostFn<'a, P> {
  graph: &'a Graph,
  partition: &'a RefCell<P>,
  branch_lengths: BTreeMap<GraphEdgeKey, f64>,
  root_key: treetime_graph::node::GraphNodeKey,
}

impl<P> GtrRateCostFn<'_, P>
where
  P: PartitionMarginalPasses + HasGtr,
{
  fn neg_log_lh(&self, sqrt_mu: f64) -> f64 {
    {
      let mut p = self.partition.borrow_mut();
      p.gtr_mut().mu = sqrt_mu * sqrt_mu;
      p.reset_node_log_likelihoods();
      if let Err(e) = marginal_backward(self.graph, &self.branch_lengths, std::slice::from_mut(&mut *p)) {
        warn!(
          "GTR rate optimization: backward pass failed at mu={:.6}: {e}",
          sqrt_mu * sqrt_mu
        );
        return f64::INFINITY;
      }
    }
    -self.partition.borrow().get_log_lh(self.root_key)
  }
}

impl<P> CostFunction for &GtrRateCostFn<'_, P>
where
  P: PartitionMarginalPasses + HasGtr,
{
  type Param = f64;
  type Output = f64;

  fn cost(&self, sqrt_mu: &Self::Param) -> Result<Self::Output, Error> {
    Ok(self.neg_log_lh(*sqrt_mu))
  }
}
