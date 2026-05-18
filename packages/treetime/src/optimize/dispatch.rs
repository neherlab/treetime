use crate::optimize::indel::estimate_indel_rate;
use crate::optimize::likelihood::evaluate_with_indels;
use crate::optimize::method_brent::{brent_inner, brent_log_inner, brent_sqrt_inner};
use crate::optimize::method_newton::{newton_inner, newton_log_inner, newton_sqrt_inner};
use crate::optimize::params::BranchOptMethod;
use crate::optimize::zero_boundary::{is_zero_branch_optimal, min_branch_length_for_indels, reconcile_zero_boundary};
use crate::representation::partition::optimization_contribution::OptimizationContribution;
use crate::representation::partition::traits::PartitionOptimizeOps;
use crate::{make_error, make_internal_report};
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

/// Unified optimization function for mixed partition types.
///
/// Main optimization loop that works with both sparse and dense partitions simultaneously.
/// For each edge, it collects contributions from all partitions and optimizes the branch
/// length using the selected method.
pub fn run_optimize_mixed<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  method: BranchOptMethod,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  P: PartitionOptimizeOps + ?Sized,
{
  let total_length = total_sequence_length(partitions);
  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot optimize branch lengths");
  }

  let indel_rate = estimate_indel_rate(graph, partitions);
  run_optimize_mixed_inner(graph, partitions, method, indel_rate, false)
}

pub fn run_optimize_mixed_with_indel_rate<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  method: BranchOptMethod,
  indel_rate: f64,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  P: PartitionOptimizeOps + ?Sized,
{
  run_optimize_mixed_inner(graph, partitions, method, indel_rate, false)
}

pub fn run_optimize_mixed_inner<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  method: BranchOptMethod,
  indel_rate: f64,
  no_indels: bool,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  P: PartitionOptimizeOps + ?Sized,
{
  let total_length = total_sequence_length(partitions);

  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot optimize branch lengths");
  }

  let one_mutation = 1.0 / total_length as f64;

  graph.get_edges().iter().try_for_each(|edge_ref| -> Result<(), Report> {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let mut branch_length = edge.branch_length().unwrap_or(0.0);

    let contributions: Vec<OptimizationContribution> = partitions
      .iter()
      .map(|partition| partition.read_arc().create_edge_contribution(edge_key))
      .collect::<Result<_, _>>()?;

    let indel_count: usize = if no_indels {
      0
    } else {
      partitions
        .iter()
        .map(|partition| partition.read_arc().edge_indel_count(edge_key))
        .sum()
    };

    // The Poisson log-likelihood derivative diverges at t=0 when k > 0, producing
    // inf/NaN in Newton's method. Use a non-zero starting point for indel-bearing edges.
    if branch_length == 0.0 && indel_count > 0 {
      branch_length = if indel_rate > 0.0 {
        (indel_count as f64 / indel_rate).max(one_mutation)
      } else {
        one_mutation
      };
    }

    // When branch length is zero and any site has non-positive likelihood at t=0
    // (mismatched certain states), evaluating ln(site_lh) or dividing by site_lh
    // produces -inf/inf. Bump to a small positive value so the evaluator operates
    // in the well-defined domain.
    if branch_length == 0.0 && !contributions.iter().all(|c| c.all_sites_valid_at_zero()) {
      branch_length = one_mutation;
    }

    // When indels are present on this edge, the Poisson derivative at t=0 is +infinity,
    // so zero branch length is never optimal. Only check the substitution-based criterion
    // when there are no indels.
    if indel_count == 0 && is_zero_branch_optimal(&contributions) {
      edge.set_branch_length(Some(0.0));
      return Ok(());
    }

    // Lower bound for Newton/Brent steps on indel-bearing edges. The Poisson derivative
    // diverges at t=0, so we must prevent the optimizer from landing exactly at zero.
    let min_branch_length = min_branch_length_for_indels(indel_count, one_mutation);

    let new_branch_length = match method {
      BranchOptMethod::Brent => brent_inner(
        branch_length,
        &contributions,
        indel_count,
        indel_rate,
        min_branch_length,
        one_mutation,
      ),
      BranchOptMethod::BrentSqrt => brent_sqrt_inner(
        branch_length,
        &contributions,
        indel_count,
        indel_rate,
        min_branch_length,
        one_mutation,
      ),
      BranchOptMethod::BrentLog => brent_log_inner(
        branch_length,
        &contributions,
        indel_count,
        indel_rate,
        min_branch_length,
        one_mutation,
      ),
      BranchOptMethod::Newton => {
        let metrics = evaluate_with_indels(&contributions, indel_count, indel_rate, branch_length);
        newton_inner(
          branch_length,
          &metrics,
          &contributions,
          indel_count,
          indel_rate,
          min_branch_length,
          one_mutation,
        )
      },
      BranchOptMethod::NewtonSqrt => {
        let metrics = evaluate_with_indels(&contributions, indel_count, indel_rate, branch_length);
        newton_sqrt_inner(
          branch_length,
          &metrics,
          &contributions,
          indel_count,
          indel_rate,
          min_branch_length,
          one_mutation,
        )
      },
      BranchOptMethod::NewtonLog => {
        // ln(t) requires t > 0; bump zero branch lengths to one_mutation
        let bl = if branch_length == 0.0 {
          one_mutation
        } else {
          branch_length
        };
        let metrics = evaluate_with_indels(&contributions, indel_count, indel_rate, bl);
        newton_log_inner(
          bl,
          &metrics,
          &contributions,
          indel_count,
          indel_rate,
          min_branch_length,
          one_mutation,
        )
      },
    }?;

    // Post-optimization boundary reconciliation. See `reconcile_zero_boundary`
    // for the full rationale. The input `branch_length` is used to size the
    // verification grid (the optimizer's output may be clamped to zero or to a
    // tiny floor like $10^{-12}$ and is unsuitable as an extent).
    let new_branch_length = reconcile_zero_boundary(
      new_branch_length,
      branch_length,
      &contributions,
      indel_count,
      indel_rate,
      one_mutation,
    )?;

    edge.set_branch_length(Some(new_branch_length));
    Ok(())
  })
}

/// Initial estimation of branch lengths for mixed partitions.
///
/// Computes per-edge substitution count over canonical (non-ambiguous,
/// non-deletion) positions via `edge_subs().len()` for both sparse and
/// dense partitions.
///
/// The denominator is the per-edge effective alignment length rather than the raw
/// sequence length, so gap-heavy edges get correctly scaled rates.
///
/// For edges with indels but no substitutions, the Poisson maximum likelihood estimate (MLE) $\hat{t} = k / \mu$
/// provides a non-zero initial estimate so Newton optimization starts from a
/// reasonable point (the indel derivative diverges at $t = 0$).
///
/// When `overwrite_valid` is false, edges that already have a finite branch
/// length are skipped, preserving calibrated input values while filling in
/// edges with missing (`None`) or invalid (`NaN`) branch lengths.
pub fn initial_guess_mixed<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  overwrite_valid: bool,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  P: PartitionOptimizeOps + ?Sized,
{
  let total_length: usize = partitions
    .iter()
    .map(|partition| partition.read_arc().sequence_length())
    .sum();

  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot compute initial guess");
  }

  let one_mutation = 1.0 / total_length as f64;
  let indel_rate = estimate_indel_rate(graph, partitions);

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.read_arc().key();
    let mut edge = edge_ref.write_arc().payload().write_arc();

    let indel_count: usize = partitions
      .iter()
      .map(|partition| partition.read_arc().edge_indel_count(edge_key))
      .sum();

    if !overwrite_valid {
      if let Some(bl) = edge.branch_length() {
        // A finite positive BL is always valid. A zero BL is valid only
        // if the edge has no indels: with indels present, the Poisson
        // derivative diverges at t=0 and estimate_indel_rate() needs
        // positive total BL to produce a nonzero rate.
        if bl.is_finite() && (bl > 0.0 || indel_count == 0) {
          continue;
        }
      }
    }

    let sub_count: usize = partitions
      .iter()
      .map(|partition| partition.read_arc().edge_subs(graph, edge_key).map(|subs| subs.len()))
      .sum::<Result<_, _>>()?;

    let effective_length: usize = partitions
      .iter()
      .map(|partition| partition.read_arc().edge_effective_length(graph, edge_key))
      .sum::<Result<_, _>>()?;

    let branch_length = if effective_length > 0 {
      let sub_estimate = sub_count as f64 / effective_length as f64;
      if sub_estimate == 0.0 && indel_count > 0 {
        if indel_rate > 0.0 {
          // Poisson MLE for indel-only branches: t = k / mu
          indel_count as f64 / indel_rate
        } else {
          // Bootstrap: rate unknown (e.g. all-zero input tree), seed with one_mutation
          // so that estimate_indel_rate produces a positive rate in run_optimize_mixed
          one_mutation
        }
      } else {
        sub_estimate
      }
    } else if indel_count > 0 {
      // All positions are gaps/ambiguous but indels are present
      one_mutation
    } else {
      one_mutation * 0.1
    };
    edge.set_branch_length(Some(branch_length));
  }

  Ok(())
}

fn total_sequence_length<P>(partitions: &[Arc<RwLock<P>>]) -> usize
where
  P: PartitionOptimizeOps + ?Sized,
{
  partitions
    .iter()
    .map(|partition| partition.read_arc().sequence_length())
    .sum()
}
