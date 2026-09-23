use crate::gtr::brent_bracketed::BrentBracketed;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, MutationCounts, infer_gtr_impl};
use crate::make_internal_report;
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalEdges, MarginalPasses, MarginalUpdate};
use argmin::core::{CostFunction, Error, Executor};
use eyre::Report;
use log::{debug, info, warn};
use ndarray::Array1;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;

pub(crate) fn refine_gtr_model<P: MarginalPasses>(
  partition: &P,
  gtr: GTR,
  update: MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>,
  iterations: usize,
  pc: f64,
  graph: &Graph,
  profile_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<(GTR, MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>), Report> {
  let MarginalUpdate {
    node_states,
    edges: MarginalEdges { backward, forward, .. },
    ..
  } = update;
  let options = InferGtrOptions {
    pc,
    ..InferGtrOptions::default()
  };

  let mut gtr = gtr;
  for i in 0..=iterations {
    let counts = partition.count_transitions(&gtr, graph, profile_lengths, &node_states, &backward, &forward)?;
    gtr = infer_gtr(&counts, &options, gtr.pi.len())?;
    debug!("GTR refinement: iteration {i}, mu = {:.6}", gtr.mu);
  }

  let update = partition.marginal_update(&gtr, graph, profile_lengths, node_states)?;
  log_final(&gtr, update.log_lh);
  Ok((gtr, update))
}

pub(crate) fn refine_gtr_model_and_rate<P: MarginalPasses>(
  partition: &P,
  gtr: GTR,
  update: MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>,
  iterations: usize,
  fixed_pi: Option<&Array1<f64>>,
  pc: f64,
  sampling_bias_correction: Option<f64>,
  graph: &Graph,
  profile_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<(GTR, MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>), Report>
where
  P::Node: Clone,
{
  let MarginalUpdate {
    node_states,
    edges: MarginalEdges { backward, forward, .. },
    ..
  } = update;
  let options = InferGtrOptions {
    fixed_pi: fixed_pi.cloned(),
    pc,
    ..InferGtrOptions::default()
  };

  let mut gtr = gtr;
  let mut nodes = node_states;
  let mut backward = backward;
  for i in 0..=iterations {
    let counts = partition.count_transitions(&gtr, graph, profile_lengths, &nodes, &backward, &forward)?;
    gtr = infer_gtr(&counts, &options, gtr.pi.len())?;
    (gtr, nodes, backward) = optimize_gtr_rate(partition, gtr, &nodes, graph, profile_lengths)?;
    debug!("GTR refinement: iteration {i}, mu = {:.6}", gtr.mu);
  }

  if let Some(correction) = sampling_bias_correction {
    gtr.mu *= correction;
    info!("Applied sampling bias correction {correction:.4}, mu = {:.6}", gtr.mu);
  }

  let update = partition.marginal_update(&gtr, graph, profile_lengths, nodes)?;
  log_final(&gtr, update.log_lh);
  Ok((gtr, update))
}

fn infer_gtr(counts: &MutationCounts, options: &InferGtrOptions, n_states: usize) -> Result<GTR, Report> {
  let result = infer_gtr_impl(counts, options)?;
  build_gtr_from_inference(n_states, &result)
}

fn build_gtr_from_inference(n_states: usize, result: &InferGtrResult) -> Result<GTR, Report> {
  GTR::builder()
    .n_states(n_states)
    .mu(result.mu)
    .W(result.W.clone())
    .pi(result.pi.clone())
    .build()
}

fn log_final(gtr: &GTR, log_lh: LogLh) {
  info!(
    "GTR refinement: final log likelihood = {:.4}, mu = {:.6}, pi = {:?}",
    log_lh.value(),
    gtr.mu,
    gtr.pi
  );
}

fn optimize_gtr_rate<P: MarginalPasses>(
  partition: &P,
  gtr: GTR,
  nodes: &BTreeMap<GraphNodeKey, P::Node>,
  graph: &Graph,
  profile_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<
  (
    GTR,
    BTreeMap<GraphNodeKey, P::Node>,
    BTreeMap<GraphEdgeKey, P::Backward>,
  ),
  Report,
>
where
  P::Node: Clone,
{
  let old_mu = gtr.mu;
  let sqrt_old_mu = old_mu.sqrt();

  let lo = 0.01 * sqrt_old_mu;
  let hi = 100.0 * sqrt_old_mu;

  let evaluate = |sqrt_mu: f64| -> (f64, Option<GtrRateCandidate<P>>) {
    let mut candidate_gtr = gtr.clone();
    candidate_gtr.mu = sqrt_mu * sqrt_mu;
    let mut candidate_nodes = nodes.clone();
    partition.reset_node_log_lh(&mut candidate_nodes);
    match partition.marginal_backward(&candidate_gtr, graph, profile_lengths, &candidate_nodes) {
      Ok(MarginalBackward {
        node_states: candidate_nodes,
        backward,
      }) => match partition.root_log_lh(graph, &candidate_nodes) {
        Ok(log_lh) => (
          -log_lh.value(),
          Some(GtrRateCandidate {
            gtr: candidate_gtr,
            nodes: candidate_nodes,
            backward,
          }),
        ),
        Err(e) => {
          warn!(
            "GTR rate optimization: root likelihood failed at mu={:.6}: {e}",
            sqrt_mu * sqrt_mu
          );
          (f64::INFINITY, None)
        },
      },
      Err(e) => {
        warn!(
          "GTR rate optimization: backward pass failed at mu={:.6}: {e}",
          sqrt_mu * sqrt_mu
        );
        (f64::INFINITY, None)
      },
    }
  };

  let neg_log_lh = |sqrt_mu: f64| evaluate(sqrt_mu).0;
  let cost_fn = GtrRateCostFn {
    neg_log_lh: &neg_log_lh,
  };

  let (cost_lo, _) = evaluate(lo);
  let (cost_mid, _) = evaluate(sqrt_old_mu);
  let (cost_hi, at_hi) = evaluate(hi);

  if cost_mid < cost_lo && cost_mid < cost_hi {
    let solver = BrentBracketed::builder().xa(lo).xb(sqrt_old_mu).xc(hi).build();
    let res = Executor::new(&cost_fn, solver)
      .configure(|cfg| cfg.max_iters(500))
      .run()
      .map_err(|e| make_internal_report!("GTR rate optimization: argmin BrentBracketed failed: {e}"))?;

    let optimal_sqrt_mu = res
      .state()
      .best_param
      .ok_or_else(|| make_internal_report!("GTR rate optimization: solver succeeded but reported no best_param"))?;

    let (_, at_opt) = evaluate(optimal_sqrt_mu);
    let GtrRateCandidate { gtr, nodes, backward } =
      at_opt.ok_or_else(|| make_internal_report!("GTR rate optimization: selected candidate failed to evaluate"))?;
    debug!(
      "GTR rate optimization: optimized mu = {:.6} (from {:.6})",
      gtr.mu, old_mu
    );
    Ok((gtr, nodes, backward))
  } else {
    let (mut restored_gtr, restored_nodes, restored_backward) =
      if let Some(GtrRateCandidate { gtr, nodes, backward }) = at_hi {
        (gtr, nodes, backward)
      } else {
        let MarginalBackward { node_states, backward } =
          partition.marginal_backward(&gtr, graph, profile_lengths, nodes)?;
        (gtr, node_states, backward)
      };
    restored_gtr.mu = old_mu;
    debug!("GTR rate optimization: skipped (no bracket), keeping mu = {old_mu:.6}");
    Ok((restored_gtr, restored_nodes, restored_backward))
  }
}

struct GtrRateCandidate<P: MarginalPasses> {
  gtr: GTR,
  nodes: BTreeMap<GraphNodeKey, P::Node>,
  backward: BTreeMap<GraphEdgeKey, P::Backward>,
}

impl<F> CostFunction for &GtrRateCostFn<'_, F>
where
  F: Fn(f64) -> f64,
{
  type Param = f64;
  type Output = f64;

  fn cost(&self, sqrt_mu: &Self::Param) -> Result<Self::Output, Error> {
    Ok((self.neg_log_lh)(*sqrt_mu))
  }
}

struct GtrRateCostFn<'a, F> {
  neg_log_lh: &'a F,
}
