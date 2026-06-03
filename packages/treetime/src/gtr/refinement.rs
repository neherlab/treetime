use crate::ancestral::marginal::{marginal_backward, update_marginal};
use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, infer_gtr_impl};
use crate::make_internal_report;
use crate::partition::traits::{HasGtr, PartitionMarginalPasses, TransitionCounting};
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::brent::BrentOpt;
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

  let cost_fn = GtrRateCostFn {
    graph,
    partition,
    partitions,
    root_key,
  };

  let cost_lo = cost_fn.neg_log_lh(lo);
  let cost_mid = cost_fn.neg_log_lh(sqrt_old_mu);
  let cost_hi = cost_fn.neg_log_lh(hi);

  if cost_mid < cost_lo && cost_mid < cost_hi {
    let solver = BrentOpt::new(lo, hi);
    let res = Executor::new(&cost_fn, solver)
      .configure(|cfg| cfg.max_iters(100))
      .run()
      .map_err(|e| make_internal_report!("GTR rate optimization: argmin BrentOpt failed: {e}"))?;

    let optimal_sqrt_mu = res
      .state()
      .best_param
      .ok_or_else(|| make_internal_report!("GTR rate optimization: solver succeeded but reported no best_param"))?;

    cost_fn.neg_log_lh(optimal_sqrt_mu);
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

struct GtrRateCostFn<'a, N: GraphNode, E: EdgeOptimizeOps, P> {
  graph: &'a Graph<N, E, ()>,
  partition: &'a Arc<RwLock<P>>,
  partitions: &'a [Arc<RwLock<P>>],
  root_key: treetime_graph::node::GraphNodeKey,
}

impl<N, E, P> GtrRateCostFn<'_, N, E, P>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + HasGtr,
{
  fn neg_log_lh(&self, sqrt_mu: f64) -> f64 {
    {
      let mut guard = self.partition.write_arc();
      guard.gtr_mut().mu = sqrt_mu * sqrt_mu;
      guard.reset_node_log_likelihoods();
    }
    match marginal_backward(self.graph, self.partitions) {
      Ok(()) => -self.partition.read_arc().get_log_lh(self.root_key),
      Err(e) => {
        warn!(
          "GTR rate optimization: backward pass failed at mu={:.6}: {e}",
          sqrt_mu * sqrt_mu
        );
        f64::INFINITY
      },
    }
  }
}

impl<N, E, P> CostFunction for &GtrRateCostFn<'_, N, E, P>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
  P: PartitionMarginalPasses<N, E> + HasGtr,
{
  type Param = f64;
  type Output = f64;

  fn cost(&self, sqrt_mu: &Self::Param) -> Result<Self::Output, Error> {
    Ok(self.neg_log_lh(*sqrt_mu))
  }
}
