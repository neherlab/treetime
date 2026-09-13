use crate::ancestral::marginal::profile_branch_lengths;
use crate::gtr::brent_bracketed::BrentBracketed;
use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::infer_gtr::common::{InferGtrOptions, InferGtrResult, infer_gtr_impl};
use crate::make_internal_report;
use crate::partition::marginal::shared::update::{
  MarginalBackward, MarginalEdges, MarginalUpdate, PartitionMarginalOps,
};
use crate::partition::traits::HasGtr;
use argmin::core::{CostFunction, Error, Executor};
use eyre::Report;
use log::{debug, info, warn};
use ndarray::Array1;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

/// Refine the GTR model by alternating inference from posterior-weighted transition counts with
/// optional substitution-rate optimization, returning the refined partition together with the marginal
/// update its final model produces.
///
/// Stable inputs (`graph`, `branch_lengths`) are borrowed; the previous marginal update is consumed and
/// a new one returned. Each rate candidate builds its own model and reconstruction on an independent
/// clone, so no candidate observes state left behind by an earlier one.
#[allow(clippy::too_many_arguments)]
pub fn refine_gtr_iterative<P>(
  graph: &Graph,
  mut partition: P,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  update: MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>,
  iterations: usize,
  fixed_pi: Option<&Array1<f64>>,
  pc: f64,
  sampling_bias_correction: Option<f64>,
  optimize_rate: bool,
) -> Result<(P, MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>), Report>
where
  P: PartitionMarginalOps + HasGtr + Clone,
  P::Node: Clone,
{
  let MarginalUpdate {
    node_states: nodes,
    edges: MarginalEdges { backward, forward, .. },
    ..
  } = update;
  let n_states = partition.gtr().pi.len();
  let options = InferGtrOptions {
    fixed_pi: fixed_pi.cloned(),
    pc,
    ..InferGtrOptions::default()
  };

  // The loop's transition counts read the message maps left by the previous step. With rate
  // optimization off (ancestral), no pass runs inside the loop, so backward/forward stay the values the
  // caller's initial update produced and only `gtr` iterates; the final update refreshes them. With
  // rate optimization on (mugration), each candidate's backward pass refreshes the backward messages.
  let mut nodes = nodes;
  let mut backward = backward;

  let counts = partition.count_transitions(graph, branch_lengths, &nodes, &backward, &forward)?;
  let result = infer_gtr_impl(&counts, &options)?;
  *partition.gtr_mut() = build_gtr_from_inference(n_states, &result)?;
  debug!("GTR refinement: initial inference, mu = {:.6}", partition.gtr().mu);

  if optimize_rate {
    (partition, nodes, backward) = optimize_gtr_rate(graph, partition, branch_lengths, &nodes)?;
    debug!(
      "GTR refinement: initial rate optimization, mu = {:.6}",
      partition.gtr().mu
    );
  }

  for i in 0..iterations {
    let counts = partition.count_transitions(graph, branch_lengths, &nodes, &backward, &forward)?;
    let result = infer_gtr_impl(&counts, &options)?;
    *partition.gtr_mut() = build_gtr_from_inference(n_states, &result)?;

    if optimize_rate {
      (partition, nodes, backward) = optimize_gtr_rate(graph, partition, branch_lengths, &nodes)?;
    }
    debug!("GTR refinement: iteration {i}, mu = {:.6}", partition.gtr().mu);
  }

  if let Some(correction) = sampling_bias_correction {
    partition.gtr_mut().mu *= correction;
    info!(
      "Applied sampling bias correction {correction:.4}, mu = {:.6}",
      partition.gtr().mu
    );
  }

  let update = partition.marginal_update(graph, &profile_branch_lengths(branch_lengths), nodes)?;

  let gtr = partition.gtr();
  info!(
    "GTR refinement: final log likelihood = {:.4}, mu = {:.6}, pi = {:?}",
    update.log_lh.value(),
    gtr.mu,
    gtr.pi
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
fn optimize_gtr_rate<P>(
  graph: &Graph,
  partition: P,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  nodes: &BTreeMap<GraphNodeKey, P::Node>,
) -> Result<(P, BTreeMap<GraphNodeKey, P::Node>, BTreeMap<GraphEdgeKey, P::Backward>), Report>
where
  P: PartitionMarginalOps + HasGtr + Clone,
  P::Node: Clone,
{
  let old_mu = partition.gtr().mu;
  let sqrt_old_mu = old_mu.sqrt();

  let lo = 0.01 * sqrt_old_mu;
  let hi = 100.0 * sqrt_old_mu;

  // Branch lengths stay fixed while only the substitution rate is optimized, so derive the profile
  // map once and route the same map into every backward pass the Brent search evaluates.
  let branch_lengths = profile_branch_lengths(branch_lengths);

  let cost_fn = GtrRateCostFn {
    graph,
    partition: &partition,
    branch_lengths: &branch_lengths,
    nodes,
  };

  let (cost_lo, _) = cost_fn.evaluate(lo);
  let (cost_mid, _) = cost_fn.evaluate(sqrt_old_mu);
  let (cost_hi, at_hi) = cost_fn.evaluate(hi);

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

    let (_, at_opt) = cost_fn.evaluate(optimal_sqrt_mu);
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
      let MarginalBackward { node_states, backward } = partition.marginal_backward(graph, &branch_lengths, nodes)?;
      (partition, node_states, backward)
    };
    restored_partition.gtr_mut().mu = old_mu;
    debug!("GTR rate optimization: skipped (no bracket), keeping mu = {old_mu:.6}");
    Ok((restored_partition, restored_nodes, restored_backward))
  }
}

/// One evaluated rate candidate: the partition at that rate together with the node states and backward
/// messages its backward pass produced.
struct GtrRateCandidate<P: PartitionMarginalOps> {
  partition: P,
  nodes: BTreeMap<GraphNodeKey, P::Node>,
  backward: BTreeMap<GraphEdgeKey, P::Backward>,
}

struct GtrRateCostFn<'a, P: PartitionMarginalOps> {
  graph: &'a Graph,
  partition: &'a P,
  branch_lengths: &'a BTreeMap<GraphEdgeKey, f64>,
  nodes: &'a BTreeMap<GraphNodeKey, P::Node>,
}

impl<P> GtrRateCostFn<'_, P>
where
  P: PartitionMarginalOps + HasGtr + Clone,
  P::Node: Clone,
{
  /// Evaluate one candidate `sqrt_mu` on an independent clone of the partition and node states: set the
  /// rate, clear the carried-over profile log likelihoods, run the backward pass, and return the
  /// negative root log likelihood together with the evaluated partition, node states, and backward
  /// messages. Clearing the log likelihoods keeps the cost the backward likelihood alone, so a forward
  /// pass's posterior log likelihood does not enter the rate search. A failed backward pass yields an
  /// infinite cost and no state, leaving the borrowed base observations intact for the next candidate.
  fn evaluate(&self, sqrt_mu: f64) -> (f64, Option<GtrRateCandidate<P>>) {
    let mut partition = self.partition.clone();
    partition.gtr_mut().mu = sqrt_mu * sqrt_mu;
    let mut nodes = self.nodes.clone();
    partition.reset_node_log_lh(&mut nodes);
    match partition.marginal_backward(self.graph, self.branch_lengths, &nodes) {
      Ok(MarginalBackward {
        node_states: nodes,
        backward,
      }) => match partition.root_log_lh(self.graph, &nodes) {
        Ok(log_lh) => (
          -log_lh.value(),
          Some(GtrRateCandidate {
            partition,
            nodes,
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
  }

  fn neg_log_lh(&self, sqrt_mu: f64) -> f64 {
    self.evaluate(sqrt_mu).0
  }
}

impl<P> CostFunction for &GtrRateCostFn<'_, P>
where
  P: PartitionMarginalOps + HasGtr + Clone,
  P::Node: Clone,
{
  type Param = f64;
  type Output = f64;

  fn cost(&self, sqrt_mu: &Self::Param) -> Result<Self::Output, Error> {
    Ok(self.neg_log_lh(*sqrt_mu))
  }
}
