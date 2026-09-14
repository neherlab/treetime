use crate::gtr::brent_bracketed::BrentBracketed;
use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, infer_gtr_impl};
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

/// Refine the GTR model at a fixed substitution rate: iterate model inference from posterior-weighted
/// transition counts, then return the refined partition together with the marginal update its final
/// model produces. This is the ancestral-reconstruction refinement, which does not optimize the rate.
///
/// The backward and forward messages from the caller's initial update stay frozen and feed every count;
/// only the model iterates. Each count still recomputes because [`MarginalPasses::count_transitions`]
/// reads the current model's `expQt`. The final update refreshes the messages under the converged model.
pub fn refine_gtr_model<P: MarginalPasses>(
  mut partition: P,
  update: MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>,
  iterations: usize,
  pc: f64,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  profile_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<(P, MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>), Report> {
  let MarginalUpdate {
    node_states,
    edges: MarginalEdges { backward, forward, .. },
    ..
  } = update;
  let options = InferGtrOptions {
    pc,
    ..InferGtrOptions::default()
  };

  for i in 0..=iterations {
    infer_and_set_gtr(
      &mut partition,
      graph,
      branch_lengths,
      &node_states,
      &backward,
      &forward,
      &options,
    )?;
    debug!("GTR refinement: iteration {i}, mu = {:.6}", partition.gtr().mu);
  }

  let update = partition.marginal_update(graph, profile_lengths, node_states)?;
  log_final(&partition, &update);
  Ok((partition, update))
}

/// Refine the GTR model and optimize the substitution rate: after each model inference, run a bracketed
/// Brent search over the rate on independent per-candidate reconstructions. This is the mugration
/// refinement.
///
/// The forward messages from the caller's initial update stay frozen and feed every count; the node
/// states and backward messages are refreshed by each rate search. `fixed_pi` pins the equilibrium
/// frequencies, and `sampling_bias_correction` scales the final rate.
pub fn refine_gtr_model_and_rate<P: MarginalPasses + Clone>(
  mut partition: P,
  update: MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>,
  iterations: usize,
  fixed_pi: Option<&Array1<f64>>,
  pc: f64,
  sampling_bias_correction: Option<f64>,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  profile_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<(P, MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>), Report>
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
  let mut nodes = node_states;
  let mut backward = backward;
  for i in 0..=iterations {
    infer_and_set_gtr(
      &mut partition,
      graph,
      branch_lengths,
      &nodes,
      &backward,
      &forward,
      &options,
    )?;
    (partition, nodes, backward) = optimize_gtr_rate(partition, &nodes, graph, profile_lengths)?;
    debug!("GTR refinement: iteration {i}, mu = {:.6}", partition.gtr().mu);
  }

  if let Some(correction) = sampling_bias_correction {
    let mut model = partition.gtr().clone();
    model.mu *= correction;
    partition.set_gtr(model);
    info!(
      "Applied sampling bias correction {correction:.4}, mu = {:.6}",
      partition.gtr().mu
    );
  }

  let update = partition.marginal_update(graph, profile_lengths, nodes)?;
  log_final(&partition, &update);
  Ok((partition, update))
}

/// Count transitions under the current model, infer a new model from the counts, and install it.
fn infer_and_set_gtr<P: MarginalPasses>(
  partition: &mut P,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  node_states: &BTreeMap<GraphNodeKey, P::Node>,
  backward: &BTreeMap<GraphEdgeKey, P::Backward>,
  forward: &BTreeMap<GraphEdgeKey, P::Forward>,
  options: &InferGtrOptions,
) -> Result<(), Report> {
  let counts = partition.count_transitions(graph, branch_lengths, node_states, backward, forward)?;
  let n_states = partition.gtr().pi.len();
  let result = infer_gtr_impl(&counts, options)?;
  partition.set_gtr(build_gtr_from_inference(n_states, &result)?);
  Ok(())
}

fn build_gtr_from_inference(n_states: usize, result: &InferGtrResult) -> Result<GTR, Report> {
  GTR::new(GTRParams {
    n_states,
    mu: result.mu,
    W: Some(result.W.clone()),
    pi: result.pi.clone(),
  })
}

fn log_final<P: MarginalPasses>(partition: &P, update: &MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>) {
  let model = partition.gtr();
  info!(
    "GTR refinement: final log likelihood = {:.4}, mu = {:.6}, pi = {:?}",
    update.log_lh.value(),
    model.mu,
    model.pi
  );
}

/// Optimize only the substitution rate `mu` by a bracketed Brent search over `sqrt(mu)`, returning the
/// partition at the selected rate together with the node states and backward messages at that rate.
///
/// Each candidate rate is evaluated on an independent clone of the partition and node states, so a
/// candidate never observes state from an earlier one. When no interior bracket is found, `mu` is
/// restored to its original value while the node states and backward messages from the last (`hi`)
/// evaluation are kept, matching the established behavior; when a bracket is found, the selected
/// candidate's evaluation supplies both the rate and the reconstruction.
fn optimize_gtr_rate<P: MarginalPasses + Clone>(
  partition: P,
  nodes: &BTreeMap<GraphNodeKey, P::Node>,
  graph: &Graph,
  profile_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<(P, BTreeMap<GraphNodeKey, P::Node>, BTreeMap<GraphEdgeKey, P::Backward>), Report>
where
  P::Node: Clone,
{
  let old_mu = partition.gtr().mu;
  let sqrt_old_mu = old_mu.sqrt();

  let lo = 0.01 * sqrt_old_mu;
  let hi = 100.0 * sqrt_old_mu;

  // Evaluate one candidate `sqrt_mu` on an independent clone of the partition and node states: set the
  // rate, clear the carried-over profile log likelihoods, run the backward pass, and return the negative
  // root log likelihood together with the evaluated partition, node states, and backward messages.
  // Clearing the log likelihoods keeps the cost the backward likelihood alone, so a forward pass's
  // posterior log likelihood does not enter the rate search. A failed backward pass yields an infinite
  // cost and no state, leaving the borrowed base observations intact for the next candidate.
  let evaluate = |sqrt_mu: f64| -> (f64, Option<GtrRateCandidate<P>>) {
    let mut candidate = partition.clone();
    let mut model = candidate.gtr().clone();
    model.mu = sqrt_mu * sqrt_mu;
    candidate.set_gtr(model);
    let mut candidate_nodes = nodes.clone();
    candidate.reset_node_log_lh(&mut candidate_nodes);
    match candidate.marginal_backward(graph, profile_lengths, &candidate_nodes) {
      Ok(MarginalBackward {
        node_states: candidate_nodes,
        backward,
      }) => match candidate.root_log_lh(graph, &candidate_nodes) {
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
      partition.gtr().mu,
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
      let MarginalBackward { node_states, backward } = partition.marginal_backward(graph, profile_lengths, nodes)?;
      (partition, node_states, backward)
    };
    let mut model = restored_partition.gtr().clone();
    model.mu = old_mu;
    restored_partition.set_gtr(model);
    debug!("GTR rate optimization: skipped (no bracket), keeping mu = {old_mu:.6}");
    Ok((restored_partition, restored_nodes, restored_backward))
  }
}

/// One evaluated rate candidate: the partition at that rate together with the node states and backward
/// messages its backward pass produced.
struct GtrRateCandidate<P: MarginalPasses> {
  partition: P,
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
