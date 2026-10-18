use crate::clock::clock_set::ClockSet;
use crate::clock::clock_state::ClockState;
use crate::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::clock::find_best_root::method_exact::optimize_exact;
use crate::clock::find_best_root::params::RootObjective;
use crate::error::input_error;
use crate::progress::LogSink;
use crate::progress_info;
use eyre::Report;
use log::debug;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

pub(crate) fn find_best_root(
  graph: &Graph,
  state: &ClockState,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  force_positive: bool,
  objective: RootObjective,
  log: &dyn LogSink,
) -> Result<FindRootResult, Report> {
  progress_info!(log, "Starting root optimization, force_positive={force_positive}");

  let root = graph.get_exactly_one_root()?;

  let root_clock_set = state.node(root.key()).clone();
  let root_acceptable = !force_positive || objective.has_positive_rate(&root_clock_set);
  let best_chisq = if root_acceptable {
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

  let edges: Vec<GraphEdgeKey> = graph.get_edges().map(|edge| edge.key()).collect();
  let candidates = edges
    .into_par_iter()
    .map(|edge| {
      optimize_exact(
        edge,
        &BranchPointCostFunction::new(state, edge, branch_lengths, objective),
      )
    })
    .collect::<Vec<_>>();
  let n_candidates = candidates.len();
  let mut rejected_negative_rate = 0;
  for res in candidates {
    if force_positive && !objective.has_positive_rate(&res.clock_set) {
      rejected_negative_rate += 1;
      continue;
    }
    if res.chisq < best_res.chisq {
      best_res = res;
    }
  }
  debug!(
    "Evaluated the best root position on {n_candidates} edges, rejected {rejected_negative_rate} with negative rate, \
     best chi-squared: {:.6e}",
    best_res.chisq
  );

  if force_positive && !objective.has_positive_rate(&best_res.clock_set) {
    return Err(input_error(
      "Clock rate is negative for all root positions. \
       The data may lack temporal signal. Please specify --clock-rate explicitly.",
    ));
  }

  progress_info!(
    log,
    "Root optimization completed. Final chi-squared: {:.6e}, split: {:.6}",
    best_res.chisq,
    best_res.split
  );

  Ok(best_res)
}

#[derive(Debug)]
pub(crate) struct FindRootResult {
  pub(crate) edge: Option<GraphEdgeKey>,

  pub(crate) split: f64,

  pub(crate) clock_set: ClockSet,

  pub(crate) chisq: f64,
}
