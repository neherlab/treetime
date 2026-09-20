use crate::clock::clock_regression::ClockVarianceParams;
use crate::clock::clock_set::ClockSet;
use crate::clock::clock_state::{ClockInputs, ClockState};
use crate::clock::find_best_root::find_best_split::{FindRootResult, find_best_split};
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RootObjective};
use crate::make_error;
use eyre::Report;
use log::{debug, info};
use rayon::prelude::*;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_utils::collections::container::get_exactly_one;

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
pub fn find_best_root(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &ClockState,
  options: &ClockVarianceParams,
  params: &BranchPointOptimizationParams,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  force_positive: bool,
  objective: RootObjective,
) -> Result<FindRootResult, Report> {
  info!("Starting root optimization with method: {params:?}, force_positive={force_positive}");

  let root = graph.get_exactly_one_root()?;
  let mut best_root_node = root;

  let root_clock_set = state.node(root.key()).clock_set.clone();
  let root_acceptable = !force_positive || has_positive_clock_rate(&root_clock_set);
  let mut best_chisq = if root_acceptable {
    objective.score(&root_clock_set)
  } else {
    f64::INFINITY
  };
  debug!(
    "Initial root chi-squared: {:.6e} (acceptable: {root_acceptable})",
    objective.score(&root_clock_set)
  );

  let mut best_res = FindRootResult {
    edge: None,
    split: 0.0,
    chisq: best_chisq,
    clock_set: root_clock_set,
  };

  let mut node_count = 0;
  let mut improvements = 0;
  let mut rejected_negative_rate = 0;
  let candidates = graph
    .get_nodes()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|node| {
      let clock_set = &state.node(node.key()).clock_set;
      let acceptable = !force_positive || has_positive_clock_rate(clock_set);
      (node, acceptable.then(|| objective.score(clock_set)))
    })
    .collect::<Vec<_>>();
  for (n, score) in candidates {
    let Some(tmp_chisq) = score else {
      rejected_negative_rate += 1;
      node_count += 1;
      continue;
    };
    if tmp_chisq < best_chisq {
      improvements += 1;
      debug!("Found better node {improvements}: chi-squared improved from {best_chisq:.6e} to {tmp_chisq:.6e}");
      best_chisq = tmp_chisq;
      best_root_node = n;
    }
    node_count += 1;
  }
  debug!(
    "Evaluated {node_count} nodes, found {improvements} improvements, \
     rejected {rejected_negative_rate} with negative rate, best chi-squared: {best_chisq:.6e}"
  );

  if !best_root_node.is_root() {
    debug!("Optimizing position on parent branch");
    let inbound = best_root_node.inbound();
    let edge = get_exactly_one(inbound).expect("Not implemented: multiple parent nodes");
    let res = find_best_split(graph, inputs, state, *edge, branch_lengths, options, params, objective)?;
    debug!(
      "Parent branch optimization result: chi-squared = {:.6e}, split = {:.6}",
      res.chisq, res.split
    );
    let split_acceptable = !force_positive || has_positive_clock_rate(&res.clock_set);
    if res.chisq < best_chisq && split_acceptable {
      debug!(
        "Parent branch optimization improved chi-squared from {:.6e} to {:.6e}",
        best_chisq, res.chisq
      );
      best_chisq = res.chisq;
      best_res = res;
    }
  }

  for (child_branch_count, e) in best_root_node.outbound().iter().enumerate() {
    debug!("Optimizing position on child branch {child_branch_count}");
    let res = find_best_split(graph, inputs, state, *e, branch_lengths, options, params, objective)?;
    debug!(
      "Child branch {} optimization result: chi-squared = {:.6e}, split = {:.6}",
      child_branch_count, res.chisq, res.split
    );
    let split_acceptable = !force_positive || has_positive_clock_rate(&res.clock_set);
    if res.chisq < best_chisq && split_acceptable {
      debug!(
        "Child branch {} optimization improved chi-squared from {:.6e} to {:.6e}",
        child_branch_count, best_chisq, res.chisq
      );
      best_chisq = res.chisq;
      best_res = res;
    }
  }

  if force_positive && !has_positive_clock_rate(&best_res.clock_set) {
    return make_error!(
      "Clock rate is negative for all root positions. \
       The data may lack temporal signal. Please specify --clock-rate explicitly."
    );
  }

  info!(
    "Root optimization completed. Final chi-squared: {:.6e}, split: {:.6}",
    best_res.chisq, best_res.split
  );

  Ok(best_res)
}

fn has_positive_clock_rate(clock_set: &ClockSet) -> bool {
  let det = clock_set.determinant();
  det > 0.0 && clock_set.clock_rate(det) > 0.0
}
