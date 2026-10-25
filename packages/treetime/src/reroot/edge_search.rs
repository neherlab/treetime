use crate::branch_lengths::branch_length_or_zero;
use crate::clock::clock_set::ClockSet;
use crate::clock::clock_state::ClockState;
use crate::error::NoPositiveRateRoot;
use crate::progress::LogSink;
use crate::progress_info;
use crate::reroot::find_root::EdgePoint;
use crate::reroot::params::RootObjective;
use eyre::Report;
use log::debug;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::split_fraction::SplitFraction;

const FLAT_CURVATURE: f64 = 1e-12;

pub(crate) fn search_edges(
  graph: &Graph,
  state: &ClockState,
  branch_lengths: &EdgeMap<Option<f64>>,
  force_positive: bool,
  objective: RootObjective,
  log: &dyn LogSink,
) -> Result<EdgeSearchResult, Report> {
  progress_info!(log, "Starting root optimization, force_positive={force_positive}");

  let root = graph.get_exactly_one_root()?;

  let root_clock_set = state.node(root.key()).clone();
  let root_acceptable = !force_positive || objective.has_positive_rate(&root_clock_set);
  let root_chisq = objective.score(&root_clock_set);
  debug!("Initial root chi-squared: {root_chisq:.6e} (acceptable: {root_acceptable})");
  let best_chisq = if root_acceptable { root_chisq } else { f64::INFINITY };

  let mut best_res = EdgeSearchResult {
    point: None,
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
    .collect::<Result<Vec<_>, Report>>()?;
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
    return Err(Report::new(NoPositiveRateRoot));
  }

  progress_info!(
    log,
    "Root optimization completed. Final chi-squared: {:.6e}, split: {:.6}",
    best_res.chisq,
    best_res.point.map_or(0.0, |point| point.split.value())
  );

  Ok(best_res)
}

#[derive(Debug)]
pub(crate) struct EdgeSearchResult {
  pub(crate) point: Option<EdgePoint>,

  pub(crate) clock_set: ClockSet,

  pub(crate) chisq: f64,
}

fn optimize_exact(edge: GraphEdgeKey, cost_fn: &BranchPointCostFunction) -> Result<EdgeSearchResult, Report> {
  let score_at = |x: SplitFraction| cost_fn.score_clock_set(&cost_fn.evaluate_clock_set(x));
  let (s0, s_half, s1) = (
    score_at(SplitFraction::PARENT_END),
    score_at(SplitFraction::MIDDLE),
    score_at(SplitFraction::CHILD_END),
  );
  let a = 2.0 * (s0 - 2.0 * s_half + s1);
  let b = s1 - s0 - a;
  let scale = s0.abs().max(s1.abs()).max(f64::MIN_POSITIVE);
  let split = SplitFraction::new(if a > FLAT_CURVATURE * scale {
    (-b / (2.0 * a)).clamp(0.0, 1.0)
  } else if s1 < s0 {
    1.0
  } else {
    0.0
  })?;
  let clock_set = cost_fn.evaluate_clock_set(split);
  let chisq = cost_fn.score_clock_set(&clock_set);
  debug!(
    "Exact root position on edge {edge}: split = {:.6}, chi-squared = {chisq:.6e}",
    split.value()
  );
  Ok(EdgeSearchResult {
    point: Some(EdgePoint { edge, split }),
    clock_set,
    chisq,
  })
}

struct BranchPointCostFunction {
  to_parent: ClockSet,
  to_child: ClockSet,
  branch_length: f64,
  objective: RootObjective,
}

impl BranchPointCostFunction {
  fn new(
    state: &ClockState,
    edge: GraphEdgeKey,
    branch_lengths: &EdgeMap<Option<f64>>,
    objective: RootObjective,
  ) -> Self {
    let edge_state = state.edge(edge);
    Self {
      to_parent: edge_state.clock_to_parent.clone(),
      to_child: edge_state.clock_to_child.clone(),
      branch_length: branch_length_or_zero(branch_lengths, edge),
      objective,
    }
  }

  fn evaluate_clock_set(&self, split: SplitFraction) -> ClockSet {
    let lengths = split.lengths(self.branch_length);
    self.to_child.shift(lengths.parent_side) + self.to_parent.shift(lengths.child_side)
  }

  fn score_clock_set(&self, clock_set: &ClockSet) -> f64 {
    self.objective.score(clock_set)
  }
}
