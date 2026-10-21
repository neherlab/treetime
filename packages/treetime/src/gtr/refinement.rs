use crate::gtr::brent_bracketed::BrentBracketed;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, MutationCounts, infer_gtr_from_counts};
use crate::make_internal_report;
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalPasses, MarginalUpdate, TransitionCounts};
use crate::progress::LogSink;
use crate::{progress_info, progress_warn};
use argmin::core::{CostFunction, Error, Executor};
use eyre::Report;
use log::debug;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;

pub(crate) fn refine_gtr<P: TransitionCounts>(
  partition: &P,
  gtr: GTR,
  input: &P::BackwardInput,
  refinement: &GtrRefinement,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<(GTR, MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>), Report>
where
  P::Node: Clone,
{
  let mut counted = partition.marginal_update(&gtr, graph, branch_lengths, input)?;
  progress_info!(
    log,
    "GTR refinement: initial log likelihood = {:.4}",
    counted.log_lh.value()
  );
  let mut gtr = gtr;
  let mut rated: Option<RatedMessages<P>> = None;
  for i in 0..=refinement.iterations {
    let backward = rated.as_ref().map_or(&counted.edges.backward, |rated| &rated.backward);
    let counts = partition.count_transitions(
      &gtr,
      graph,
      branch_lengths,
      &counted.node_states,
      backward,
      &counted.edges.forward,
    )?;
    gtr = infer_gtr(&counts, &refinement.inference, gtr.pi().len(), log)?;
    if rated.is_none() {
      counted = partition.marginal_update(&gtr, graph, branch_lengths, P::backward_input(&counted.node_states))?;
      debug!(
        "GTR refinement: log likelihood under the first fitted model = {:.4}",
        counted.log_lh.value()
      );
    }
    let nodes = rated.as_ref().map_or(&counted.node_states, |rated| &rated.nodes);
    let (rated_gtr, nodes, backward) = optimize_gtr_rate(
      partition,
      gtr,
      nodes,
      graph,
      branch_lengths,
      RateRangeEdge::KeepStart,
      log,
    )?;
    gtr = rated_gtr;
    rated = Some(RatedMessages { nodes, backward });
    debug!("GTR refinement: iteration {i}, mu = {:.6}", gtr.mu());
  }

  if let Some(correction) = refinement.sampling_bias_correction {
    gtr.set_mu(gtr.mu() * correction)?;
    progress_info!(
      log,
      "Applied sampling bias correction {correction:.4}, mu = {:.6}",
      gtr.mu()
    );
  }

  let nodes = rated.as_ref().map_or(&counted.node_states, |rated| &rated.nodes);
  let update = partition.marginal_update(&gtr, graph, branch_lengths, P::backward_input(nodes))?;
  log_final(&gtr, update.log_lh, log);
  Ok((gtr, update))
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct GtrRefinement {
  pub(crate) iterations: usize,
  pub(crate) inference: InferGtrOptions,
  pub(crate) sampling_bias_correction: Option<f64>,
}

struct RatedMessages<P: MarginalPasses> {
  nodes: BTreeMap<GraphNodeKey, P::Node>,
  backward: BTreeMap<GraphEdgeKey, P::Backward>,
}

fn infer_gtr(
  counts: &MutationCounts,
  options: &InferGtrOptions,
  n_states: usize,
  log: &dyn LogSink,
) -> Result<GTR, Report> {
  let result = infer_gtr_from_counts(counts, options, log)?;
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

fn log_final(gtr: &GTR, log_lh: LogLh, log: &dyn LogSink) {
  progress_info!(
    log,
    "GTR refinement: final log likelihood = {:.4}, mu = {:.6}, pi = {:?}",
    log_lh.value(),
    gtr.mu(),
    gtr.pi()
  );
}

pub(crate) fn optimize_gtr_rate<P: MarginalPasses>(
  partition: &P,
  gtr: GTR,
  nodes: &BTreeMap<GraphNodeKey, P::Node>,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  range_edge: RateRangeEdge,
  log: &dyn LogSink,
) -> Result<RateStepResult<P>, Report>
where
  P::Node: Clone,
{
  let old_mu = gtr.mu();
  let sqrt_old_mu = old_mu.sqrt();

  let lo = 0.01 * sqrt_old_mu;
  let hi = 100.0 * sqrt_old_mu;

  let evaluate = |sqrt_mu: f64| -> (f64, Option<GtrRateCandidate<P>>) {
    let mut candidate_gtr = gtr.clone();
    if let Err(e) = candidate_gtr.set_mu(sqrt_mu * sqrt_mu) {
      progress_warn!(log, "GTR rate optimization: invalid rate at sqrt(mu)={sqrt_mu:.6}: {e}");
      return (f64::INFINITY, None);
    }
    let input = P::fresh_backward_input(nodes);
    match partition.marginal_backward(&candidate_gtr, graph, branch_lengths, &input) {
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
          progress_warn!(
            log,
            "GTR rate optimization: root likelihood failed at mu={:.6}: {e}",
            sqrt_mu * sqrt_mu
          );
          (f64::INFINITY, None)
        },
      },
      Err(e) => {
        progress_warn!(
          log,
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

  let (cost_lo, at_lo) = evaluate(lo);
  let (cost_mid, at_mid) = evaluate(sqrt_old_mu);
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
      gtr.mu(),
      old_mu
    );
    Ok((gtr, nodes, backward))
  } else {
    match range_edge {
      RateRangeEdge::KeepStart => keep_start_rate(partition, gtr, at_mid, nodes, graph, branch_lengths),
      RateRangeEdge::LowerCostEnd => {
        let (cost_end, at_end) = if cost_lo <= cost_hi {
          (cost_lo, at_lo)
        } else {
          (cost_hi, at_hi)
        };
        let GtrRateCandidate { gtr, nodes, backward } = at_end.ok_or_else(|| {
          make_internal_report!("GTR rate optimization: both ends of the rate range failed to evaluate")
        })?;
        debug!(
          "GTR rate optimization: no bracket, taking mu = {:.6e} (from {old_mu:.6}) at the range end with cost {cost_end:.6}",
          gtr.mu()
        );
        Ok((gtr, nodes, backward))
      },
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RateRangeEdge {
  KeepStart,
  LowerCostEnd,
}

fn keep_start_rate<P: MarginalPasses>(
  partition: &P,
  gtr: GTR,
  at_start: Option<GtrRateCandidate<P>>,
  nodes: &BTreeMap<GraphNodeKey, P::Node>,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<RateStepResult<P>, Report> {
  let old_mu = gtr.mu();
  let (mut restored_gtr, restored_nodes, restored_backward) =
    if let Some(GtrRateCandidate { gtr, nodes, backward }) = at_start {
      (gtr, nodes, backward)
    } else {
      let MarginalBackward { node_states, backward } =
        partition.marginal_backward(&gtr, graph, branch_lengths, P::backward_input(nodes))?;
      (gtr, node_states, backward)
    };
  restored_gtr.set_mu(old_mu)?;
  debug!("GTR rate optimization: skipped (no bracket), keeping mu = {old_mu:.6}");
  Ok((restored_gtr, restored_nodes, restored_backward))
}

type RateStepResult<P> = (
  GTR,
  BTreeMap<GraphNodeKey, <P as MarginalPasses>::Node>,
  BTreeMap<GraphEdgeKey, <P as MarginalPasses>::Backward>,
);

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
