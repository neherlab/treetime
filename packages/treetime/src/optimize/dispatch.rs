#![allow(
  clippy::as_conversions,
  reason = "counts and indices to f64 for normalization and coordinate math"
)]

use crate::optimize::branch_length::{is_valid_branch_length_value, validate_branch_length_value};
use crate::optimize::indel::estimate_indel_rate;
use crate::optimize::likelihood::evaluate_with_indels;
use crate::optimize::method_brent::{brent_inner, brent_log_inner, brent_sqrt_inner};
use crate::optimize::method_newton::{newton_inner, newton_log_inner, newton_sqrt_inner};
use crate::optimize::params::BranchOptMethod;
use crate::optimize::zero_boundary::{is_zero_branch_optimal, min_branch_length_for_indels, reconcile_zero_boundary};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::{make_error, make_internal_report, make_report};
use eyre::{Report, WrapErr};
use rayon::prelude::*;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

/// Unified optimization function for mixed partition types.
///
/// Main optimization loop that works with both sparse and dense partitions simultaneously. For each
/// edge, it reads the pre-gathered per-partition contributions and indel count and optimizes the branch
/// length using the selected method.
pub fn run_optimize_mixed(
  graph: &Graph,
  total_length: usize,
  contributions: &BTreeMap<GraphEdgeKey, Vec<OptimizationContribution>>,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  method: BranchOptMethod,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot optimize branch lengths");
  }

  let indel_rate = estimate_indel_rate(graph, indel_counts, branch_lengths);
  run_optimize_mixed_inner(
    graph,
    total_length,
    contributions,
    indel_counts,
    method,
    indel_rate,
    false,
    branch_lengths,
  )?;
  Ok(())
}

#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub fn run_optimize_mixed_with_indel_rate(
  graph: &Graph,
  total_length: usize,
  contributions: &BTreeMap<GraphEdgeKey, Vec<OptimizationContribution>>,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  method: BranchOptMethod,
  indel_rate: f64,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  run_optimize_mixed_inner(
    graph,
    total_length,
    contributions,
    indel_counts,
    method,
    indel_rate,
    false,
    branch_lengths,
  )?;
  Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn run_optimize_mixed_inner(
  graph: &Graph,
  total_length: usize,
  contributions: &BTreeMap<GraphEdgeKey, Vec<OptimizationContribution>>,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  method: BranchOptMethod,
  indel_rate: f64,
  no_indels: bool,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot optimize branch lengths");
  }

  let one_mutation = 1.0 / total_length as f64;

  // For a bifurcating root under reversible models, only the sum of root-edge branch
  // lengths is identifiable (Pulley Principle). Capture the pre-optimization ratio of the
  // two root edges; after the per-edge loop, redistribute the optimized total in this ratio.
  // See kb/issues/M-optimize-root-bifurcating-independent-vs-joint.md
  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.key();
    let branch_length = branch_lengths[&edge_key]
      .ok_or_else(|| make_report!("Cannot optimize edge {edge_key} with a missing branch length"))?;
    validate_branch_length_value(branch_length).wrap_err_with(|| format!("Cannot optimize edge {edge_key}"))?;
  }

  let root_state = BifurcatingRootState::capture(graph, branch_lengths)?;

  // Each per-edge optimum depends only on that edge's input length and its partition
  // contributions, so the parallel map is order free; collect the results and write them into the
  // shared branch-length map serially afterward (a `BTreeMap` cannot be mutated from parallel
  // workers).
  let branch_lengths_in = &*branch_lengths;
  let optimized: Vec<(GraphEdgeKey, f64)> = graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| -> Result<(GraphEdgeKey, f64), Report> {
      let edge_key = edge_ref.key();
      let mut branch_length = branch_lengths_in[&edge_key]
        .ok_or_else(|| make_internal_report!("Validated edge {edge_key} lost its branch length"))?;

      let contributions = &contributions[&edge_key];

      let indel_count: usize = if no_indels { 0 } else { indel_counts[&edge_key] };

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
      if indel_count == 0 && is_zero_branch_optimal(contributions) {
        return Ok((edge_key, 0.0));
      }

      // Lower bound for Newton/Brent steps on indel-bearing edges. The Poisson derivative
      // diverges at t=0, so we must prevent the optimizer from landing exactly at zero.
      let min_branch_length = min_branch_length_for_indels(indel_count, one_mutation);

      let new_branch_length = match method {
        BranchOptMethod::Brent => brent_inner(
          branch_length,
          contributions,
          indel_count,
          indel_rate,
          min_branch_length,
          one_mutation,
        ),
        BranchOptMethod::BrentSqrt => brent_sqrt_inner(
          branch_length,
          contributions,
          indel_count,
          indel_rate,
          min_branch_length,
          one_mutation,
        ),
        BranchOptMethod::BrentLog => brent_log_inner(
          branch_length,
          contributions,
          indel_count,
          indel_rate,
          min_branch_length,
          one_mutation,
        ),
        BranchOptMethod::Newton => {
          let metrics = evaluate_with_indels(contributions, indel_count, indel_rate, branch_length)?;
          newton_inner(
            branch_length,
            &metrics,
            contributions,
            indel_count,
            indel_rate,
            min_branch_length,
            one_mutation,
          )
        },
        BranchOptMethod::NewtonSqrt => {
          let metrics = evaluate_with_indels(contributions, indel_count, indel_rate, branch_length)?;
          newton_sqrt_inner(
            branch_length,
            &metrics,
            contributions,
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
          let metrics = evaluate_with_indels(contributions, indel_count, indel_rate, bl)?;
          newton_log_inner(
            bl,
            &metrics,
            contributions,
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
        contributions,
        indel_count,
        indel_rate,
        one_mutation,
      )?;

      Ok((edge_key, new_branch_length))
    })
    .collect::<Result<Vec<_>, Report>>()?;

  for (edge_key, new_branch_length) in optimized {
    branch_lengths.insert(edge_key, Some(new_branch_length));
  }

  if let Some(state) = root_state {
    state.restore(branch_lengths);
  }

  Ok(())
}

struct BifurcatingRootState {
  edge0: GraphEdgeKey,
  edge1: GraphEdgeKey,
  ratio: f64,
}

impl BifurcatingRootState {
  fn capture(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Result<Option<Self>, Report> {
    let root = graph.get_exactly_one_root()?;
    let children = graph.children_of(root).collect::<Vec<_>>();
    if children.len() == 2 {
      let edge0 = children[0].1.key();
      let edge1 = children[1].1.key();
      let bl0 = branch_lengths[&edge0].unwrap_or(0.0);
      let bl1 = branch_lengths[&edge1].unwrap_or(0.0);
      let total = bl0 + bl1;
      let ratio = if total > 0.0 { bl0 / total } else { 0.5 };
      Ok(Some(Self { edge0, edge1, ratio }))
    } else {
      Ok(None)
    }
  }

  fn restore(self, branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>) {
    let Self { edge0, edge1, ratio } = self;
    let total = branch_lengths[&edge0].unwrap_or(0.0) + branch_lengths[&edge1].unwrap_or(0.0);
    branch_lengths.insert(edge0, Some(total * ratio));
    branch_lengths.insert(edge1, Some(total * (1.0 - ratio)));
  }
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
///
/// When `no_indels` is true, indel counts and rates do not affect either
/// branch validity or the estimated branch length.
#[allow(clippy::too_many_arguments)]
pub fn initial_guess_mixed(
  graph: &Graph,
  total_length: usize,
  indel_counts: &BTreeMap<GraphEdgeKey, usize>,
  sub_counts: &BTreeMap<GraphEdgeKey, usize>,
  effective_lengths: &BTreeMap<GraphEdgeKey, usize>,
  overwrite_valid: bool,
  no_indels: bool,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  if total_length == 0 {
    return make_error!("Total sequence length across all partitions is zero; cannot compute initial guess");
  }

  let one_mutation = 1.0 / total_length as f64;
  let indel_rate = if no_indels {
    0.0
  } else {
    estimate_indel_rate(graph, indel_counts, branch_lengths)
  };

  for edge_ref in graph.get_edges() {
    let edge_key = edge_ref.key();

    let indel_count: usize = if no_indels { 0 } else { indel_counts[&edge_key] };

    if !overwrite_valid {
      if let Some(bl) = branch_lengths.get(&edge_key).copied().flatten() {
        // A finite positive BL is always valid. A zero BL is valid only
        // if the edge has no indels: with indels present, the Poisson
        // derivative diverges at t=0 and estimate_indel_rate() needs
        // positive total BL to produce a nonzero rate.
        if is_valid_branch_length_value(bl) && (bl > 0.0 || indel_count == 0) {
          continue;
        }
      }
    }

    let sub_count: usize = sub_counts[&edge_key];

    let effective_length: usize = effective_lengths[&edge_key];

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
    branch_lengths.insert(edge_key, Some(branch_length));
  }

  Ok(())
}
