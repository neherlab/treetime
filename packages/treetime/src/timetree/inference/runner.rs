use crate::branch_lengths::one_mutation;
use crate::clock::clock_model::ClockModel;
use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::mean_field::CoalescentPrior;
use crate::optimize::gather::EdgeContributions;
use crate::optimize::indel::IndelPolicy;
use crate::partition::set::PartitionSet;
use crate::progress::LogSink;
use crate::timetree::branch_model::BranchModel;
use crate::timetree::inference::backward_pass::propagate_distributions_backward;
use crate::timetree::inference::bad_branches::derive_bad_branches;
use crate::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
use crate::timetree::inference::forward_pass::{ForwardPassInputs, propagate_distributions_forward};
use crate::timetree::inference::result::{
  BranchLikelihood, NodePosterior, NodeTimes, PosteriorSummaries, TimeInference, log_marginal_likelihood,
};
use crate::{progress_info, progress_warn};
use eyre::Report;
use log::debug;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_grid::{GridPointLimitExceeded, MaxGridPoints};
use treetime_primitives::LogLh;
use treetime_utils::fmt::float::float_to_significant_digits;

pub(crate) const GRID_POINTS: usize = 300;

const _: () = assert!(
  GRID_POINTS < MaxGridPoints::MIN,
  "every grid point limit must admit the grids of GRID_POINTS points that every run creates"
);

pub(crate) const EPS: f64 = 5e-4;

pub(crate) const TRIM_EPS: f64 = 5e-4;

pub(crate) fn run_timetree(
  inputs: &TimeInferenceInputs<'_>,
  coalescent: Option<&CoalescentPrior>,
  log: &dyn LogSink,
) -> Result<TimeInference, Report> {
  infer_times(inputs, coalescent, log)
    .map_err(|report| explain_grid_point_limit(report, inputs.clock_model.clock_rate(), inputs.clock_rate_fixed))
}

fn infer_times(
  inputs: &TimeInferenceInputs<'_>,
  coalescent: Option<&CoalescentPrior>,
  log: &dyn LogSink,
) -> Result<TimeInference, Report> {
  let TimeInferenceInputs {
    graph,
    date_constraints,
    leaf_bad_branches,
    gammas,
    branch_model,
    branch_lengths,
    names,
    clock_model,
    indel_policy,
    max_grid_points,
    posterior_summaries,
    ..
  } = *inputs;
  progress_info!(log, "# Running timetree inference");

  progress_info!(log, "## Using clock model");
  let clock_rate = clock_model.clock_rate();
  progress_info!(log, "**Clock rate:** {clock_rate:.6e}");

  let bad_branches = derive_bad_branches(graph, date_constraints, leaf_bad_branches)?;

  let branches = match branch_model {
    BranchModel::Input => {
      progress_info!(log, "## Creating branch distributions from input lengths");
      create_branch_distributions_input_mode(graph, branch_lengths, gammas, clock_rate)
    },
    BranchModel::Marginal(partitions) => {
      progress_info!(log, "## Computing branch distributions from partitions");
      compute_branch_distributions_marginal_mode(
        graph,
        partitions,
        branch_lengths,
        gammas,
        clock_rate,
        indel_policy,
        max_grid_points,
        log,
      )?
    },
  };

  progress_info!(log, "## Propagating distributions backward");
  let backward = propagate_distributions_backward(
    graph,
    date_constraints,
    coalescent,
    &bad_branches,
    &branches,
    max_grid_points,
  )?;

  let log_lh_time = log_marginal_likelihood(&backward, graph.root_key()?)?.map(LogLh::new);
  if let Some(log_lh) = log_lh_time {
    progress_info!(log, "**Log marginal likelihood of node times:** {:.4}", log_lh.value());
  }

  progress_info!(log, "## Propagating distributions forward");
  let forward_inputs = ForwardPassInputs {
    constraints: date_constraints,
    names,
    branches: &branches,
    max_points: max_grid_points,
  };
  let posterior = propagate_distributions_forward(
    graph,
    &forward_inputs,
    backward,
    |key, view| NodePosterior::summarize(key, view, posterior_summaries, bad_branches[&key]),
    log,
  )?;

  progress_info!(log, "# Timetree inference completed");
  Ok(TimeInference {
    log_lh_time,
    bad_branches,
    branches,
    posterior,
    posterior_summaries,
  })
}

pub(crate) struct TimeInferenceInputs<'a> {
  pub graph: &'a Graph,
  pub date_constraints: &'a DateConstraints,
  pub leaf_bad_branches: &'a BTreeMap<GraphNodeKey, bool>,
  pub gammas: &'a BTreeMap<GraphEdgeKey, f64>,
  pub branch_model: &'a BranchModel,
  pub branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  pub names: &'a BTreeMap<GraphNodeKey, Option<String>>,
  pub clock_model: &'a ClockModel,
  pub clock_rate_fixed: bool,
  pub indel_policy: IndelPolicy,
  pub max_grid_points: MaxGridPoints,
  pub posterior_summaries: PosteriorSummaries,
}

pub(super) fn explain_grid_point_limit(report: Report, clock_rate: f64, clock_rate_fixed: bool) -> Report {
  let Some(exceeded) = report
    .chain()
    .find_map(|cause| cause.downcast_ref::<GridPointLimitExceeded>())
    .copied()
  else {
    return report;
  };
  let limit = exceeded.limit;
  let needs = match exceeded.required_points() {
    Some(points) => format!("a grid of {points} points, more than"),
    None => "a grid with more points than".to_owned(),
  };
  let rate = float_to_significant_digits(clock_rate, 3);
  let cause = if clock_rate_fixed {
    format!(
      "The fixed --clock-rate {rate} may be far from the rate the data supports, which makes branch-time \
       distributions very narrow. Check --clock-rate"
    )
  } else {
    format!(
      "The estimated clock rate {rate} may be far from the true rate, which makes branch-time distributions very \
       narrow; errors in the input dates and weak temporal signal cause such estimates. Check the input dates, set \
       --clock-rate"
    )
  };
  report.wrap_err(format!(
    "Time inference needs {needs} the limit of {limit} (--max-grid-points). {cause}, or raise --max-grid-points \
     when enough memory is available."
  ))
}

pub(crate) const CLOCK_BRANCH_LENGTH_DAMPING: f64 = 0.5;

pub(crate) const CLOCK_BRANCH_LENGTH_UNDAMPED: f64 = 1.0;

pub(crate) fn blended_clock_branch_lengths(
  graph: &Graph,
  clock_rate: f64,
  damping: f64,
  previous_lengths: &BTreeMap<GraphEdgeKey, f64>,
  node_times: &NodeTimes,
  gammas: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> BTreeMap<GraphEdgeKey, f64> {
  let committed: Vec<(GraphEdgeKey, f64, bool)> = graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .filter_map(|edge_ref| {
      let key = edge_ref.key();
      let (Some(parent_time), Some(child_time)) = (node_times[&edge_ref.source()], node_times[&edge_ref.target()])
      else {
        return None;
      };

      let duration = child_time - parent_time;
      let fresh = clock_rate * gammas[&key] * duration.max(0.0);
      let value = match previous_lengths.get(&key) {
        Some(previous) => (1.0 - damping) * previous + damping * fresh,
        None => fresh,
      };

      Some((key, value, duration < 0.0))
    })
    .collect();

  let mut clock_branch_lengths: BTreeMap<GraphEdgeKey, f64> = graph
    .get_edges()
    .filter_map(|edge| {
      let key = edge.key();
      previous_lengths.get(&key).map(|length| (key, *length))
    })
    .collect();
  let mut inverted = 0_usize;
  for (key, value, is_inverted) in committed {
    clock_branch_lengths.insert(key, value);
    inverted += usize::from(is_inverted);
  }

  if inverted > 0 {
    progress_warn!(
      log,
      "Timetree: {inverted} branch(es) run backwards in time, i.e. the child is dated before its \
       parent. Their clock branch lengths were committed as zero. This is expected only where an \
       observed leaf date or an exact internal-node date conflicts with the fitted clock, since the \
       forward pass clamps the other internal nodes to their parent but leaves exact dates as given."
    );
  }
  clock_branch_lengths
}

#[expect(
  clippy::too_many_arguments,
  reason = "the tree, reconstruction, per-edge rate multipliers, clock rate, indel policy and grid limit are independent inputs of one time-inference pass"
)]
fn compute_branch_distributions_marginal_mode(
  graph: &Graph,
  partitions: &PartitionSet,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  gammas: &BTreeMap<GraphEdgeKey, f64>,
  clock_rate: f64,
  indel_policy: IndelPolicy,
  max_grid_points: MaxGridPoints,
  log: &dyn LogSink,
) -> Result<BTreeMap<GraphEdgeKey, BranchLikelihood>, Report> {
  let total_sites = partitions.total_length();
  let one_mutation = one_mutation(total_sites);

  let indels = indel_policy.edge_indels(graph, partitions, branch_lengths);

  let partition_count = partitions.len();
  progress_info!(
    log,
    "Computing branch distributions from {partition_count} partition(s) with {total_sites} total sites"
  );
  debug!("One mutation = {one_mutation:.6e} substitutions/site");
  let indel_rate = indels.rate();
  debug!("Indel rate = {indel_rate:.6e} indels/(site*time)");

  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| -> Result<(GraphEdgeKey, BranchLikelihood), Report> {
      let edge_key = edge_ref.key();
      let branch_length = branch_lengths[&edge_key].unwrap_or(one_mutation);
      let gamma = gammas[&edge_key];

      debug!("Edge {edge_key:?}: input branch_length = {branch_length:.6e}, gamma = {gamma:.4}");

      let distribution = partitions.with_edge(edge_key, |contributions| {
        compute_branch_length_distribution(
          contributions,
          indels.term(edge_key),
          branch_length,
          one_mutation,
          GRID_POINTS,
          clock_rate * gamma,
          max_grid_points,
        )
      })?;

      let time_length = distribution.likely_time()?;
      if let Some(likely_time) = time_length {
        debug!("Edge {edge_key:?}: distribution peak at time = {likely_time:.6e}");
      }

      let branch = BranchLikelihood {
        distribution: Some(distribution),
        time_length,
      };
      Ok((edge_key, branch))
    })
    .collect()
}

pub(super) fn create_branch_distributions_input_mode(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  gammas: &BTreeMap<GraphEdgeKey, f64>,
  clock_rate: f64,
) -> BTreeMap<GraphEdgeKey, BranchLikelihood> {
  graph
    .get_edges()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|edge_ref| {
      let key = edge_ref.key();
      // TODO: this is wrong. The branch length distribution should be a gamma distribution with branch_length/one_mutation
      // as the shape parameter. n_mut = branch_length/one_mutation --> P(dt) = (mu*dt)^n_mut * exp(-mu*dt) / n_mut!
      // v0 input mode uses a Poisson approximation for short branches and a Gaussian one for longer branches,
      // see kb/algo/unimplemented.md, section "Branch Length Interpolator (Input Mode)".
      let time_length = branch_lengths[&key].map(|branch_length| branch_length / (clock_rate * gammas[&key]));
      let branch = BranchLikelihood {
        distribution: time_length.map(|time_length| Arc::new(Distribution::point(time_length, 0.0))),
        time_length,
      };
      (key, branch)
    })
    .collect()
}

pub(crate) fn timetree_branch_lengths(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  clock_branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> BTreeMap<GraphEdgeKey, f64> {
  graph
    .get_edges()
    .map(|edge| {
      let key = edge.key();
      let branch_length = clock_branch_lengths
        .get(&key)
        .copied()
        .or_else(|| branch_lengths[&key])
        .unwrap_or(0.0);
      (key, branch_length)
    })
    .collect()
}
