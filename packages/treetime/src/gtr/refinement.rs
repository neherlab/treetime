use crate::gtr::brent_bracketed::BrentBracketed;
use crate::gtr::gtr::{GTR, GTRParams};
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

/// Refine the GTR model at a fixed substitution rate: iterate model inference from posterior-weighted
/// transition counts, then return the refined model together with the marginal update it produces. This
/// is the ancestral-reconstruction refinement, which does not optimize the rate.
///
/// The partition is an immutable source; the model flows through as a value. The backward and forward
/// messages from the caller's initial update stay frozen and feed every count; only the model iterates.
/// Each count still recomputes because [`MarginalPasses::count_transitions`] reads the current model's
/// `expQt`. The final update refreshes the messages under the converged model.
pub fn refine_gtr_model<P: MarginalPasses>(
  partition: &P,
  gtr: GTR,
  update: MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>,
  iterations: usize,
  pc: f64,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
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
    let counts = partition.count_transitions(&gtr, graph, branch_lengths, &node_states, &backward, &forward)?;
    gtr = infer_gtr(&counts, &options, gtr.pi.len())?;
    debug!("GTR refinement: iteration {i}, mu = {:.6}", gtr.mu);
  }

  let update = partition.marginal_update(&gtr, graph, profile_lengths, node_states)?;
  log_final(&gtr, update.log_lh);
  Ok((gtr, update))
}

/// Refine the GTR model and optimize the substitution rate: after each model inference, run a bracketed
/// Brent search over the rate on independent per-candidate reconstructions. This is the mugration
/// refinement.
///
/// The forward messages from the caller's initial update stay frozen and feed every count; the node
/// states and backward messages are refreshed by each rate search. `fixed_pi` pins the equilibrium
/// frequencies, and `sampling_bias_correction` scales the final rate.
pub fn refine_gtr_model_and_rate<P: MarginalPasses>(
  partition: &P,
  gtr: GTR,
  update: MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>,
  iterations: usize,
  fixed_pi: Option<&Array1<f64>>,
  pc: f64,
  sampling_bias_correction: Option<f64>,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
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

  // The count reads the frozen forward messages and the node states and backward messages from the
  // previous rate search; the final update refreshes them under the converged model.
  let mut gtr = gtr;
  let mut nodes = node_states;
  let mut backward = backward;
  for i in 0..=iterations {
    let counts = partition.count_transitions(&gtr, graph, branch_lengths, &nodes, &backward, &forward)?;
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

/// Infer a GTR model from transition counts: the representation-independent estimation core.
fn infer_gtr(counts: &MutationCounts, options: &InferGtrOptions, n_states: usize) -> Result<GTR, Report> {
  let result = infer_gtr_impl(counts, options)?;
  build_gtr_from_inference(n_states, &result)
}

fn build_gtr_from_inference(n_states: usize, result: &InferGtrResult) -> Result<GTR, Report> {
  GTR::new(GTRParams {
    n_states,
    mu: result.mu,
    W: Some(result.W.clone()),
    pi: result.pi.clone(),
  })
}

fn log_final(gtr: &GTR, log_lh: LogLh) {
  info!(
    "GTR refinement: final log likelihood = {:.4}, mu = {:.6}, pi = {:?}",
    log_lh.value(),
    gtr.mu,
    gtr.pi
  );
}

/// Optimize only the substitution rate `mu` by a bracketed Brent search over `sqrt(mu)`, returning the
/// model at the selected rate together with the node states and backward messages at that rate.
///
/// Each candidate rate is evaluated on an independent clone of the model and node states over the shared
/// immutable partition, so a candidate never observes state from an earlier one. When no interior
/// bracket is found, `mu` is restored to its original value while the node states and backward messages
/// from the last (`hi`) evaluation are kept, matching the established behavior; when a bracket is found,
/// the selected candidate's evaluation supplies both the rate and the reconstruction.
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

  // Evaluate one candidate `sqrt_mu` on an independent clone of the model and node states: set the rate,
  // clear the carried-over profile log likelihoods, run the backward pass, and return the negative root
  // log likelihood together with the evaluated model, node states, and backward messages. Clearing the
  // log likelihoods keeps the cost the backward likelihood alone, so a forward pass's posterior log
  // likelihood does not enter the rate search. A failed backward pass yields an infinite cost and no
  // state, leaving the borrowed base observations intact for the next candidate.
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
    let GtrRateCandidate { gtr, nodes, backward } =
      at_opt.ok_or_else(|| make_internal_report!("GTR rate optimization: selected candidate failed to evaluate"))?;
    debug!(
      "GTR rate optimization: optimized mu = {:.6} (from {:.6})",
      gtr.mu, old_mu
    );
    Ok((gtr, nodes, backward))
  } else {
    // No interior bracket: keep the node states and backward messages from the last (`hi`) evaluation
    // but restore the rate, exactly as before. A failed `hi` evaluation leaves the input observations
    // untouched, so fall back to the input node states with the rate restored.
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

/// One evaluated rate candidate: the model at that rate together with the node states and backward
/// messages its backward pass produced.
struct GtrRateCandidate<P: MarginalPasses> {
  gtr: GTR,
  nodes: BTreeMap<GraphNodeKey, P::Node>,
  backward: BTreeMap<GraphEdgeKey, P::Backward>,
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
