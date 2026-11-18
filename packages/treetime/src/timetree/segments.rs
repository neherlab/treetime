use crate::cancel::Cancel;
use crate::clock::clock_model::ClockModel;
use crate::clock::date_constraints::load_date_constraints;
use crate::error::OperationError;
use crate::optimize::indel::IndelPolicy;
use crate::optimize::params::{BranchOptMethod, InitialGuessMode, OptimizeLoopSettings, TopologyOps};
use crate::optimize::partition_rates::PartitionRates;
use crate::optimize::pipeline::{OptimizeInput, OptimizeParams, OptimizedGraph, optimize_graph, optimize_results};
use crate::partition::result::{AlignmentResult, PartitionResult};
use crate::partition::split::SequenceInput;
use crate::progress::{LogSink, StageSink};
use crate::progress_info;
use crate::timetree::branch_model::BranchModel;
use crate::timetree::confidence::{NodeConfidenceInterval, extract_confidence_intervals};
use crate::timetree::inference::backward_pass::ParentMessages;
use crate::timetree::inference::bad_branches::bad_leaves;
use crate::timetree::inference::result::{PosteriorSummaries, date_branch_lengths};
use crate::timetree::inference::runner::{
  CLOCK_BRANCH_LENGTH_DAMPING, TimeInferenceInputs, blended_clock_branch_lengths, run_timetree_with,
  timetree_branch_lengths,
};
use eyre::WrapErr;
use std::collections::BTreeSet;
use std::sync::Arc;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node_names::NodeNames;
use treetime_graph::scope::EdgeScopes;
use treetime_grid::MaxGridPoints;
use treetime_ops::FftPlans;
use treetime_primitives::LogLh;
use treetime_primitives::date::DateConstraint;
use treetime_utils::make_report;

const OPTIMIZE_MAX_ITER: usize = 10;

const OPTIMIZE_TOLERANCE: f64 = 0.1;

pub struct SegmentTimetreeParams {
  pub dense: Option<bool>,
  pub indel_policy: IndelPolicy,
  pub partition_rates: PartitionRates,
  pub clock_rate: f64,
  pub confidence: bool,
  pub max_grid_points: MaxGridPoints,
  pub cross_parent_messages: bool,
  pub refinement_rounds: usize,
}

pub struct SegmentTimetreeInput {
  pub graph: Graph,
  pub names: NodeNames,
  pub sequences: SequenceInput,
  pub dates: NodeMap<DateConstraint>,
  pub branch_lengths: EdgeMap<Option<f64>>,
  pub scopes: Arc<EdgeScopes>,
}

pub struct SegmentTimetreeOutput {
  pub graph: Graph,
  pub names: NodeNames,
  pub branch_lengths: EdgeMap<Option<f64>>,
  pub node_dates: NodeMap<Option<f64>>,
  pub date_branch_lengths: EdgeMap<Option<f64>>,
  pub confidence_intervals: Option<Vec<NodeConfidenceInterval>>,
  pub log_lh_time: Option<LogLh>,
  pub alignments: Vec<AlignmentResult>,
  pub partitions: Vec<PartitionResult>,
}

pub fn run_segment_timetree(
  params: &SegmentTimetreeParams,
  input: SegmentTimetreeInput,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<SegmentTimetreeOutput, OperationError> {
  progress_info!(log, "# TreeTime Timetree Estimation on a graph of segment trees");
  let clock_model = ClockModel::at_fixed_rate(params.clock_rate).map_err(OperationError::InvalidParams)?;
  let SegmentTimetreeInput {
    graph,
    names,
    sequences,
    dates,
    branch_lengths,
    scopes,
  } = input;
  let date_constraints = load_date_constraints(&dates, &graph, log)
    .wrap_err("Failed to load date constraints")
    .map_err(OperationError::InvalidInput)?;
  let optimize_input = OptimizeInput {
    graph,
    names,
    sequences,
    branch_lengths,
    scopes: Some(scopes),
  };
  let OptimizedGraph {
    graph,
    names,
    branch_lengths,
    partitions,
    alignments,
    model_names,
    indel_rates,
  } = optimize_graph(&segment_optimize_params(params), optimize_input, cancel, stages, log)?;

  cancel.check().map_err(OperationError::classify)?;
  stages.report("Time inference", 0.6, "");
  let leaf_bad_branches = bad_leaves(&graph, &date_constraints, &BTreeSet::new());
  if leaf_bad_branches.values().all(|&bad| bad) {
    return Err(OperationError::InvalidInput(make_report!(
      "None of the tips of the segment graph has a usable date"
    )));
  }
  let gammas: EdgeMap<f64> = graph.get_edges().map(|edge| (edge.key(), 1.0)).collect();
  let fft_plans = FftPlans::new();
  let parent_messages = if params.cross_parent_messages {
    ParentMessages::AllEdges
  } else {
    ParentMessages::ChildrenOnly
  };
  let mut branch_model = BranchModel::Marginal(partitions);
  let mut clock_lengths = EdgeMap::new();
  let mut round = 0;
  let inference = loop {
    let last = round == params.refinement_rounds;
    let time_inputs = TimeInferenceInputs {
      graph: &graph,
      date_constraints: &date_constraints,
      leaf_bad_branches: &leaf_bad_branches,
      gammas: &gammas,
      branch_model: &branch_model,
      branch_lengths: &branch_lengths,
      names: &names,
      clock_model: &clock_model,
      clock_rate_fixed: true,
      indel_policy: params.indel_policy,
      max_grid_points: params.max_grid_points,
      fft_plans: &fft_plans,
      posterior_summaries: PosteriorSummaries {
        regions: last && params.confidence,
        lineage_survival: false,
      },
    };
    let inference = run_timetree_with(&time_inputs, None, parent_messages, log)
      .wrap_err("Timetree inference on the segment graph failed")
      .map_err(OperationError::classify)?;
    if last {
      break inference;
    }
    round += 1;
    cancel.check().map_err(OperationError::classify)?;
    progress_info!(
      log,
      "## Refinement round {round}: marginal reconstruction at the clock branch lengths"
    );
    clock_lengths = blended_clock_branch_lengths(
      &graph,
      params.clock_rate,
      CLOCK_BRANCH_LENGTH_DAMPING,
      &clock_lengths,
      &inference.node_times(),
      &gammas,
      log,
    );
    branch_model = branch_model
      .marginal_update(
        &graph,
        &timetree_branch_lengths(&graph, &branch_lengths, &clock_lengths),
      )
      .map_err(OperationError::classify)?;
  };
  let BranchModel::Marginal(partitions) = branch_model else {
    return Err(OperationError::classify(make_report!(
      "The timetree of segment trees lost its marginal reconstruction"
    )));
  };
  let (alignments, partition_results) = optimize_results(&graph, alignments, &partitions, &model_names, &indel_rates)
    .map_err(OperationError::InferenceFailed)?;

  let node_dates = inference.node_times();
  let confidence_intervals = params
    .confidence
    .then(|| extract_confidence_intervals(&graph, &inference, &NodeMap::new(), &names))
    .transpose()
    .map_err(OperationError::InferenceFailed)?;
  Ok(SegmentTimetreeOutput {
    date_branch_lengths: date_branch_lengths(&graph, &node_dates),
    log_lh_time: inference.log_lh_time,
    graph,
    names,
    branch_lengths,
    node_dates,
    confidence_intervals,
    alignments,
    partitions: partition_results,
  })
}

fn segment_optimize_params(params: &SegmentTimetreeParams) -> OptimizeParams {
  OptimizeParams {
    dense: params.dense,
    initial_guess: InitialGuessMode::Auto,
    reroot_spec: None,
    loop_settings: OptimizeLoopSettings {
      max_iter: OPTIMIZE_MAX_ITER,
      dp: OPTIMIZE_TOLERANCE,
      damping: 0.0,
      opt_method: BranchOptMethod::Brent,
      indel_policy: params.indel_policy,
      topology_ops: TopologyOps::NONE,
      partition_rates: params.partition_rates,
    },
  }
}
