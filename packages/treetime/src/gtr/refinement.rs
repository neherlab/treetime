use crate::gtr::brent_bracketed::BrentBracketed;
use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, MutationCounts, infer_gtr_impl};
use crate::make_internal_report;
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalEdges, MarginalUpdate};
use argmin::core::{CostFunction, Error, Executor};
use eyre::Report;
use log::{debug, info, warn};
use ndarray::Array1;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;

/// Refine the GTR model by alternating inference from posterior-weighted transition counts with
/// optional substitution-rate optimization, returning the refined partition together with the marginal
/// update its final model produces.
///
/// The backend-invariant estimation core (transition-count inference, model construction, and the
/// bracketed rate search) is shared here; the representation supplies its passes and its model access
/// as closures. Each rate candidate builds its own model and reconstruction on an independent clone, so
/// no candidate observes state left behind by an earlier one.
///
/// Closures: `gtr`/`set_gtr` read and replace the partition's substitution model; `count_transitions`,
/// `run_update`, `run_backward`, and `root_log_lh` run the representation's marginal passes and read the
/// root likelihood; `reset` zeros the per-node profile log likelihoods before a rate-search backward
/// pass.
#[allow(clippy::too_many_arguments)]
pub fn refine_gtr_iterative<P, N, B, F, E>(
  mut partition: P,
  update: MarginalUpdate<N, B, F, E>,
  iterations: usize,
  fixed_pi: Option<&Array1<f64>>,
  pc: f64,
  sampling_bias_correction: Option<f64>,
  optimize_rate: bool,
  gtr: impl Fn(&P) -> &GTR,
  set_gtr: impl Fn(&mut P, GTR),
  count_transitions: impl Fn(
    &P,
    &BTreeMap<GraphNodeKey, N>,
    &BTreeMap<GraphEdgeKey, B>,
    &BTreeMap<GraphEdgeKey, F>,
  ) -> Result<MutationCounts, Report>,
  run_update: impl Fn(&P, BTreeMap<GraphNodeKey, N>) -> Result<MarginalUpdate<N, B, F, E>, Report>,
  run_backward: impl Fn(&P, &BTreeMap<GraphNodeKey, N>) -> Result<MarginalBackward<N, B>, Report>,
  root_log_lh: impl Fn(&P, &BTreeMap<GraphNodeKey, N>) -> Result<LogLh, Report>,
  reset: impl Fn(&P, &mut BTreeMap<GraphNodeKey, N>),
) -> Result<(P, MarginalUpdate<N, B, F, E>), Report>
where
  P: Clone,
  N: Clone,
{
  let MarginalUpdate {
    node_states: nodes,
    edges: MarginalEdges { backward, forward, .. },
    ..
  } = update;
  let n_states = gtr(&partition).pi.len();
  let options = InferGtrOptions {
    fixed_pi: fixed_pi.cloned(),
    pc,
    ..InferGtrOptions::default()
  };

  // The loop's transition counts read the message maps left by the previous step. With rate
  // optimization off (ancestral), no pass runs inside the loop, so backward/forward stay the values the
  // caller's initial update produced and only the model iterates; the final update refreshes them. With
  // rate optimization on (mugration), each candidate's backward pass refreshes the backward messages.
  let mut nodes = nodes;
  let mut backward = backward;

  let counts = count_transitions(&partition, &nodes, &backward, &forward)?;
  let result = infer_gtr_impl(&counts, &options)?;
  set_gtr(&mut partition, build_gtr_from_inference(n_states, &result)?);
  debug!("GTR refinement: initial inference, mu = {:.6}", gtr(&partition).mu);

  if optimize_rate {
    (partition, nodes, backward) =
      optimize_gtr_rate(partition, &nodes, &gtr, &set_gtr, &run_backward, &root_log_lh, &reset)?;
    debug!(
      "GTR refinement: initial rate optimization, mu = {:.6}",
      gtr(&partition).mu
    );
  }

  for i in 0..iterations {
    let counts = count_transitions(&partition, &nodes, &backward, &forward)?;
    let result = infer_gtr_impl(&counts, &options)?;
    set_gtr(&mut partition, build_gtr_from_inference(n_states, &result)?);

    if optimize_rate {
      (partition, nodes, backward) =
        optimize_gtr_rate(partition, &nodes, &gtr, &set_gtr, &run_backward, &root_log_lh, &reset)?;
    }
    debug!("GTR refinement: iteration {i}, mu = {:.6}", gtr(&partition).mu);
  }

  if let Some(correction) = sampling_bias_correction {
    let mut model = gtr(&partition).clone();
    model.mu *= correction;
    set_gtr(&mut partition, model);
    info!(
      "Applied sampling bias correction {correction:.4}, mu = {:.6}",
      gtr(&partition).mu
    );
  }

  let update = run_update(&partition, nodes)?;

  let model = gtr(&partition);
  info!(
    "GTR refinement: final log likelihood = {:.4}, mu = {:.6}, pi = {:?}",
    update.log_lh.value(),
    model.mu,
    model.pi
  );

  Ok((partition, update))
}

fn build_gtr_from_inference(n_states: usize, result: &InferGtrResult) -> Result<GTR, Report> {
  GTR::new(GTRParams {
    n_states,
    mu: result.mu,
    W: Some(result.W.clone()),
    pi: result.pi.clone(),
  })
}

/// Optimize only the substitution rate `mu` by a bracketed Brent search over `sqrt(mu)`, returning the
/// partition at the selected rate together with the node states and backward messages at that rate.
///
/// Each candidate rate is evaluated on an independent clone of the partition and node states, so a
/// candidate never observes state from an earlier one. When no interior bracket is found, `mu` is
/// restored to its original value while the node states and backward messages from the last (`hi`)
/// evaluation are kept, matching the established behavior; when a bracket is found, the selected
/// candidate's evaluation supplies both the rate and the reconstruction.
#[allow(clippy::too_many_arguments)]
fn optimize_gtr_rate<P, N, B>(
  partition: P,
  nodes: &BTreeMap<GraphNodeKey, N>,
  gtr: impl Fn(&P) -> &GTR,
  set_gtr: impl Fn(&mut P, GTR),
  run_backward: impl Fn(&P, &BTreeMap<GraphNodeKey, N>) -> Result<MarginalBackward<N, B>, Report>,
  root_log_lh: impl Fn(&P, &BTreeMap<GraphNodeKey, N>) -> Result<LogLh, Report>,
  reset: impl Fn(&P, &mut BTreeMap<GraphNodeKey, N>),
) -> Result<(P, BTreeMap<GraphNodeKey, N>, BTreeMap<GraphEdgeKey, B>), Report>
where
  P: Clone,
  N: Clone,
{
  let old_mu = gtr(&partition).mu;
  let sqrt_old_mu = old_mu.sqrt();

  let lo = 0.01 * sqrt_old_mu;
  let hi = 100.0 * sqrt_old_mu;

  // Evaluate one candidate `sqrt_mu` on an independent clone of the partition and node states: set the
  // rate, clear the carried-over profile log likelihoods, run the backward pass, and return the negative
  // root log likelihood together with the evaluated partition, node states, and backward messages.
  // Clearing the log likelihoods keeps the cost the backward likelihood alone, so a forward pass's
  // posterior log likelihood does not enter the rate search. A failed backward pass yields an infinite
  // cost and no state, leaving the borrowed base observations intact for the next candidate.
  let evaluate = |sqrt_mu: f64| -> (f64, Option<GtrRateCandidate<P, N, B>>) {
    let mut candidate = partition.clone();
    let mut model = gtr(&candidate).clone();
    model.mu = sqrt_mu * sqrt_mu;
    set_gtr(&mut candidate, model);
    let mut candidate_nodes = nodes.clone();
    reset(&candidate, &mut candidate_nodes);
    match run_backward(&candidate, &candidate_nodes) {
      Ok(MarginalBackward {
        node_states: candidate_nodes,
        backward,
      }) => match root_log_lh(&candidate, &candidate_nodes) {
        Ok(log_lh) => (
          -log_lh.value(),
          Some(GtrRateCandidate {
            partition: candidate,
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

    let (_, at_opt) = evaluate(optimal_sqrt_mu);
    let GtrRateCandidate {
      partition,
      nodes,
      backward,
    } = at_opt.ok_or_else(|| make_internal_report!("GTR rate optimization: selected candidate failed to evaluate"))?;
    debug!(
      "GTR rate optimization: optimized mu = {:.6} (from {:.6})",
      gtr(&partition).mu,
      old_mu
    );
    Ok((partition, nodes, backward))
  } else {
    // No interior bracket: keep the node states and backward messages from the last (`hi`) evaluation
    // but restore the rate, exactly as before. A failed `hi` evaluation leaves the input observations
    // untouched, so fall back to the input node states with the rate restored.
    let (mut restored_partition, restored_nodes, restored_backward) = if let Some(GtrRateCandidate {
      partition,
      nodes,
      backward,
    }) = at_hi
    {
      (partition, nodes, backward)
    } else {
      let MarginalBackward { node_states, backward } = run_backward(&partition, nodes)?;
      (partition, node_states, backward)
    };
    let mut model = gtr(&restored_partition).clone();
    model.mu = old_mu;
    set_gtr(&mut restored_partition, model);
    debug!("GTR rate optimization: skipped (no bracket), keeping mu = {old_mu:.6}");
    Ok((restored_partition, restored_nodes, restored_backward))
  }
}

/// One evaluated rate candidate: the partition at that rate together with the node states and backward
/// messages its backward pass produced.
struct GtrRateCandidate<P, N, B> {
  partition: P,
  nodes: BTreeMap<GraphNodeKey, N>,
  backward: BTreeMap<GraphEdgeKey, B>,
}

/// Adapter presenting the rate search's negative-log-likelihood closure to argmin's cost interface.
struct GtrRateCostFn<'a, F> {
  neg_log_lh: &'a F,
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
