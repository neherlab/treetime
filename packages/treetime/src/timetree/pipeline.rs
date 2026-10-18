use crate::alphabet::alphabet::Alphabet;
use crate::branch_lengths::branch_length_or_zero;
use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockFit, ClockTree, ClockVarianceParams};
use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
use crate::clock::divergence::root_to_node_divergences;
use crate::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::clock::gaussian::gaussian_clock;
use crate::clock::reroot::RerootParams;
use crate::clock::rtt::{ClockDateSource, ClockRegressionResult, clock_fit_regression_results};
use crate::coalescent::mean_field::CoalescentPrior;
use crate::coalescent::population_size::effective_population_size;
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::optimize::params::BranchLengthMode;
use crate::partition::create::{Representation, build_marginal_partition};
use crate::progress::{LogSink, StageSink};
use crate::seq::alignment::NodeSeqInput;
use crate::seq::mutation::{Mutation, MutationTrack, SequenceMutations, edge_state_change_counts};
use crate::seq::sink::SeqSink;
use crate::timetree::branch_model::BranchModel;
use crate::timetree::coalescent::CoalescentOutput;
use crate::timetree::coalescent_timescale::{
  CoalescentMode, CoalescentSetup, CoalescentTimescale, build_coalescent_output,
};
use crate::timetree::confidence::{
  NodeConfidenceInterval, RateSusceptibility, compute_rate_susceptibility, determine_rate_std,
  extract_confidence_intervals,
};
use crate::timetree::convergence::optimizer::TraceSink;
use crate::timetree::inference::result::given_times;
use crate::timetree::inference::runner::timetree_branch_lengths;
use crate::timetree::optimization::reroot::{DatedClockFit, fit_clock_to_dates};
use crate::timetree::params::{
  FinalSequences, TimeMarginalMode, TimetreeContext, TimetreeParams, build_covariation_clock_params,
  compute_effective_time_marginal,
};
use crate::timetree::pre_loop::{PreLoopInputs, PreLoopState, run_pre_loop};
use crate::timetree::refinement_loop::run_refinement_loop;
use crate::timetree::round::{RoundInputs, RoundState, final_marginal_round, infer_final_times, run_initial_round};
use crate::{progress_info, progress_warn};
use eyre::{Report, WrapErr};
use log::debug;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;
use treetime_primitives::date::DateConstraint;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::make_report;

pub fn run(
  params: &TimetreeParams,
  input: TimetreeInput,
  trace_sink: Option<&mut dyn TraceSink>,
  seq_sink: Option<&mut dyn SeqSink>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<TimetreeOutput, OperationError> {
  progress_info!(log, "# TreeTime Timetree Estimation");
  validate_params(params, seq_sink.is_some())?;
  let TimetreeInput {
    graph,
    alphabet,
    sequences,
    dates,
    branch_lengths,
    names,
  } = input;
  let context = prepare_inputs(params, &graph, sequences.as_ref(), dates.as_ref(), log)?;

  cancel.check().map_err(OperationError::classify)?;
  stages.report("Clock regression", 0.1, "");
  let (graph, branch_lengths, clock_fit) =
    estimate_initial_clock(params, &context, graph, branch_lengths, &names, log).map_err(OperationError::classify)?;
  let has_alignment = sequences.is_some();
  let branch_model = initialize_branch_model(params, &graph, &branch_lengths, alphabet, sequences, &names, log)?;
  let pre_loop_inputs = PreLoopInputs {
    params,
    context: &context,
    names: &names,
    has_alignment,
  };
  let pre_loop_state = PreLoopState::new(graph, branch_lengths, branch_model, clock_fit);
  let pre_loop =
    run_pre_loop(&pre_loop_inputs, pre_loop_state, cancel, stages, log).map_err(OperationError::classify)?;

  let initial = run_initial_round(params, &context, &names, pre_loop, log).map_err(OperationError::classify)?;
  let round_inputs = RoundInputs {
    params,
    context: &context,
    leaf_bad_branches: &initial.leaf_bad_branches,
    outliers: &initial.outliers,
  };
  let coalescent = &initial.coalescent;
  let (state, timescale) = run_refinement_loop(
    &round_inputs,
    coalescent,
    initial.timescale,
    initial.state,
    trace_sink,
    cancel,
    stages,
    log,
  )
  .map_err(OperationError::classify)?;
  report_coalescent_size(params, coalescent.mode, &timescale, log);

  cancel.check().map_err(OperationError::classify)?;
  stages.report("Postprocessing", 0.85, "");
  progress_info!(log, "### TreeTime: postprocessing");
  let final_times =
    refine_final_times(&round_inputs, coalescent, &timescale, state, log).map_err(OperationError::classify)?;
  let results =
    gather_results(params, &context, coalescent, &timescale, final_times).map_err(OperationError::classify)?;
  assemble_output(params, &context, seq_sink, initial.outliers, dates, results, log)
}

fn assemble_output(
  params: &TimetreeParams,
  context: &TimetreeContext,
  seq_sink: Option<&mut dyn SeqSink>,
  outliers: BTreeSet<GraphNodeKey>,
  dates: Option<BTreeMap<GraphNodeKey, DateConstraint>>,
  results: FinalResults,
  log: &dyn LogSink,
) -> Result<TimetreeOutput, OperationError> {
  let RoundState {
    graph,
    names,
    branch_model,
    branch_lengths,
    clock_model,
    clock_branch_lengths,
    time_inference,
    ..
  } = results.state;
  let gtr = branch_model.gtr().cloned();
  let model_name = branch_model.is_marginal().then_some(params.model);
  let sequences = reconstruct_final_sequences(
    params,
    context,
    seq_sink,
    &graph,
    branch_model,
    &timetree_branch_lengths(&graph, &branch_lengths, &clock_branch_lengths),
    log,
  )?;
  let node_dates = time_inference.node_times();
  let date_branch_lengths = date_branch_lengths(&graph, &node_dates);

  Ok(TimetreeOutput {
    clock_model,
    clock_regression: results.clock_regression,
    confidence_intervals: results.confidence_intervals,
    dates,
    gtr,
    model_name,
    coalescent: results.coalescent_output,
    branch_lengths,
    date_branch_lengths,
    node_dates,
    bad_branches: time_inference.bad_branches,
    divergences: results.divergences,
    outliers,
    sequences,
    names,
    graph,
  })
}

pub struct TimetreeInput {
  pub graph: Graph,
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  pub alphabet: Alphabet,
  pub sequences: Option<BTreeMap<GraphNodeKey, NodeSeqInput>>,
  pub dates: Option<BTreeMap<GraphNodeKey, DateConstraint>>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

pub struct TimetreeOutput {
  pub graph: Graph,
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  pub node_dates: BTreeMap<GraphNodeKey, Option<f64>>,
  pub divergences: BTreeMap<GraphNodeKey, f64>,
  pub outliers: BTreeSet<GraphNodeKey>,
  pub bad_branches: BTreeMap<GraphNodeKey, bool>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub date_branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub clock_model: ClockModel,
  pub clock_regression: Vec<ClockRegressionResult>,
  pub confidence_intervals: Option<Vec<NodeConfidenceInterval>>,
  pub gtr: Option<GTR>,
  pub model_name: Option<GtrModelName>,
  pub coalescent: Option<CoalescentOutput>,
  pub dates: Option<BTreeMap<GraphNodeKey, DateConstraint>>,
  pub sequences: Option<TimetreeSequences>,
}

fn validate_params(params: &TimetreeParams, has_seq_sink: bool) -> Result<(), OperationError> {
  if params.n_branches_posterior.is_some() {
    return Err(OperationError::InvalidParams(make_report!(
      "--n-branches-posterior is not yet implemented"
    )));
  }
  if has_seq_sink && params.branch_length_mode == BranchLengthMode::Input {
    return Err(OperationError::InvalidParams(make_report!(
      "Reconstructed sequence output requires ancestral reconstruction; \
       incompatible with --branch-length-mode=input"
    )));
  }
  if has_seq_sink && !params.sequence_outputs_requested {
    return Err(OperationError::InvalidParams(make_report!(
      "A sequence sink was passed, but the parameters request no sequence outputs"
    )));
  }
  Ok(())
}

fn prepare_inputs(
  params: &TimetreeParams,
  graph: &Graph,
  sequences: Option<&BTreeMap<GraphNodeKey, NodeSeqInput>>,
  dates: Option<&BTreeMap<GraphNodeKey, DateConstraint>>,
  log: &dyn LogSink,
) -> Result<TimetreeContext, OperationError> {
  debug!(
    "Branch length mode: {:?}, Keep root: {}",
    params.branch_length_mode, params.keep_root
  );

  let time_marginal = compute_effective_time_marginal(params.time_marginal, params.confidence, log);

  let covariation_clock_params = build_covariation_clock_params(
    params.covariation,
    params.sequence_length,
    params.tip_slack,
    sequences,
    log,
  )
  .map_err(OperationError::InvalidParams)?;

  let date_constraints = if let Some(dates) = dates {
    load_date_constraints(dates, graph, log)
      .wrap_err("Failed to load date constraints")
      .map_err(OperationError::InvalidInput)?
  } else {
    DateConstraints::default()
  };

  let final_sequences = if !params.sequence_outputs_requested || params.branch_length_mode != BranchLengthMode::Marginal
  {
    FinalSequences::Skip
  } else if time_marginal.runs_final_round() {
    FinalSequences::FromLastRound
  } else {
    FinalSequences::AfterMarginalUpdate
  };
  Ok(TimetreeContext {
    final_sequences,
    time_marginal,
    date_constraints,
    covariation_clock_params: covariation_clock_params.unwrap_or_default(),
    branch_params: BranchPointOptimizationParams::default(),
  })
}

fn estimate_initial_clock(
  params: &TimetreeParams,
  context: &TimetreeContext,
  graph: Graph,
  branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<(Graph, BTreeMap<GraphEdgeKey, Option<f64>>, ClockFit), Report> {
  let reroot_params = RerootParams::new(params.reroot_spec.clone(), !params.allow_negative_rate);
  let fit = DatedClockFit {
    outliers: &BTreeSet::new(),
    clock_params: &ClockVarianceParams::default(),
    clock_rate: params.clock_rate,
    keep_root: params.keep_root,
    branch_params: &context.branch_params,
    reroot_params: &reroot_params,
    failure: "Failed to infer clock model",
  };
  let (
    ClockTree {
      graph, branch_lengths, ..
    },
    clock_reroot,
  ) = fit_clock_to_dates(graph, branch_lengths, &context.date_constraints, &fit, names, log)?;
  Ok((graph, branch_lengths, clock_reroot.into_clock_fit()?))
}

fn initialize_branch_model(
  params: &TimetreeParams,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  alphabet: Alphabet,
  aln: Option<BTreeMap<GraphNodeKey, NodeSeqInput>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<BranchModel, OperationError> {
  match params.branch_length_mode {
    BranchLengthMode::Input => {
      progress_info!(log, "Branch length mode: Input - using tree branch lengths");
      Ok(BranchModel::Input)
    },
    BranchLengthMode::Marginal => {
      progress_info!(
        log,
        "Branch length mode: Marginal - initializing partitions from alignment"
      );
      let node_inputs = aln
        .ok_or_else(|| OperationError::InvalidInput(make_report!("Alignment required for marginal reconstruction")))?;
      let reconstruction = build_marginal_partition(
        Representation::resolve(params.dense),
        params.model,
        graph,
        alphabet,
        node_inputs,
        &branch_lengths_or_zero(branch_lengths),
        log,
      )
      .map_err(OperationError::classify)?;
      Ok(BranchModel::Marginal(reconstruction))
    },
  }
}

fn report_coalescent_size(
  params: &TimetreeParams,
  mode: CoalescentMode,
  timescale: &CoalescentTimescale,
  log: &dyn LogSink,
) {
  if mode.output_mode().is_none() {
    return;
  }
  let tc_values = timescale.schedule.values();
  progress_info!(
    log,
    "Coalescent effective population size (gen_per_year={:.4}, {} segment(s)):",
    params.gen_per_year,
    tc_values.len()
  );
  for (i, &tc) in tc_values.iter().enumerate() {
    progress_info!(
      log,
      "  segment {i}: Tc = {tc:.6e}  N_e = {:.6e}",
      effective_population_size(tc, params.gen_per_year)
    );
  }
}

/// Rate uncertainty of the Gaussian timetree at the final rate, or `None` with a warning when the sequence length is
/// unknown or the approximation degenerates
fn gaussian_rate_std(inputs: &RoundInputs<'_>, state: &RoundState, log: &dyn LogSink) -> Result<Option<f64>, Report> {
  let sequence_length = match state.branch_model.sequence_length() {
    0 => inputs.params.sequence_length.unwrap_or(0),
    length => length,
  };
  if sequence_length == 0 {
    progress_warn!(
      log,
      "Cannot estimate the clock rate uncertainty without a sequence length. \
       Specify --clock-std-dev or --sequence-length. Proceeding without rate uncertainty."
    );
    return Ok(None);
  }
  let tip_dates = given_times(&state.graph, &inputs.context.date_constraints)?;
  match gaussian_clock(
    &state.graph,
    &state.branch_lengths,
    sequence_length.to_f64_exact(),
    &tip_dates,
    inputs.leaf_bad_branches,
    state.clock_model.clock_rate(),
  ) {
    Ok(clock) => {
      progress_info!(
        log,
        "Gaussian clock: rate std {:.4e}, root time std {:.4e}",
        clock.rate_std,
        clock.root_time_std
      );
      Ok(Some(clock.rate_std))
    },
    Err(report) => {
      progress_warn!(log, "Proceeding without rate uncertainty: {report}");
      Ok(None)
    },
  }
}

fn refine_final_times(
  inputs: &RoundInputs<'_>,
  coalescent: &CoalescentSetup,
  timescale: &CoalescentTimescale,
  state: RoundState,
  log: &dyn LogSink,
) -> Result<FinalTimes, Report> {
  let params = inputs.params;
  let rate_std = if params.confidence {
    determine_rate_std(params.clock_std_dev, || gaussian_rate_std(inputs, &state, log))?
  } else {
    None
  };

  let final_model = CoalescentPrior::new(&coalescent.lineage_expectation, &timescale.distribution)?;
  let final_prior = coalescent.prior_wanted().then_some(&final_model);

  let (central, rate_susceptibility_dates) = if let Some(rate_std) = rate_std {
    progress_info!(log, "### Rate susceptibility analysis (rate_std={rate_std:.6e})");
    let RateSusceptibility { dates, central } =
      compute_rate_susceptibility(&state.time_inference_inputs(inputs), final_prior, rate_std, log)
        .wrap_err("Rate susceptibility analysis failed")?;
    (Some(central), dates)
  } else {
    (None, BTreeMap::new())
  };

  let state = if inputs.context.time_marginal.runs_final_round() {
    progress_info!(log, "### Final round: marginal reconstruction for confidence intervals");
    let time_inference = match central {
      Some(central) => central,
      None => infer_final_times(inputs, final_prior, &state, log)?,
    };
    final_marginal_round(time_inference, state, log)?
  } else {
    match central {
      Some(time_inference) => RoundState {
        time_inference,
        ..state
      },
      None => state,
    }
  };

  Ok(FinalTimes {
    state,
    rate_std,
    rate_susceptibility_dates,
  })
}

fn gather_results(
  params: &TimetreeParams,
  context: &TimetreeContext,
  coalescent: &CoalescentSetup,
  timescale: &CoalescentTimescale,
  final_times: FinalTimes,
) -> Result<FinalResults, Report> {
  let FinalTimes {
    state,
    rate_std,
    rate_susceptibility_dates,
  } = final_times;

  let confidence_intervals = (matches!(
    context.time_marginal,
    TimeMarginalMode::OnlyFinal | TimeMarginalMode::Always
  ) || rate_std.is_some())
  .then(|| {
    extract_confidence_intervals(
      &state.graph,
      &state.time_inference.posterior,
      &rate_susceptibility_dates,
      &state.names,
    )
  })
  .transpose()?;

  let coalescent_output = build_coalescent_output(
    coalescent.mode,
    timescale,
    params.gen_per_year,
    &coalescent.skyline_params,
  )?;

  let divergences = root_to_node_divergences(&state.graph, |edge_key| {
    branch_length_or_zero(&state.branch_lengths, edge_key)
  })?;
  let names = assign_node_names(state.names, &state.graph)?.names;

  let clock_regression = clock_fit_regression_results(&state.clock_model, &state.clock_points, &names, |key| {
    if context.date_constraints.date_constraint(key).is_some() {
      ClockDateSource::Input
    } else {
      ClockDateSource::Inferred
    }
  });

  Ok(FinalResults {
    state: RoundState { names, ..state },
    confidence_intervals,
    coalescent_output,
    divergences,
    clock_regression,
  })
}

struct FinalTimes {
  state: RoundState,
  rate_std: Option<f64>,
  rate_susceptibility_dates: BTreeMap<GraphNodeKey, [f64; 3]>,
}

struct FinalResults {
  state: RoundState,
  confidence_intervals: Option<Vec<NodeConfidenceInterval>>,
  coalescent_output: Option<CoalescentOutput>,
  divergences: BTreeMap<GraphNodeKey, f64>,
  clock_regression: Vec<ClockRegressionResult>,
}

fn reconstruct_final_sequences(
  params: &TimetreeParams,
  context: &TimetreeContext,
  mut seq_sink: Option<&mut dyn SeqSink>,
  graph: &Graph,
  branch_model: BranchModel,
  final_branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<Option<TimetreeSequences>, OperationError> {
  let BranchModel::Marginal(reconstruction) = branch_model else {
    if params.include_leaves || params.impute_missing_data {
      progress_warn!(
        log,
        "Ignoring tip-state flags (--include-leaves / --impute-missing-data / --reconstruct-tip-states): \
         no ancestral reconstruction was performed under --branch-length-mode=input"
      );
    }
    return Ok(None);
  };
  let reconstruction = match context.final_sequences {
    FinalSequences::Skip => return Ok(None),
    FinalSequences::FromLastRound => reconstruction,
    FinalSequences::AfterMarginalUpdate => {
      reconstruction
        .marginal_update(graph, final_branch_lengths)
        .map_err(OperationError::classify)?
        .0
    },
  };
  if let Some(sink) = seq_sink.as_deref_mut() {
    sink.on_topology(graph).map_err(OperationError::SinkFailed)?;
  }
  let SequenceMutations {
    root_sequence,
    edge_mutations,
  } = reconstruction.stream_sequences(
    graph,
    params.impute_missing_data,
    &MutationTrack::Nucleotide,
    params.include_leaves,
    params.report_ambiguous,
    seq_sink,
  )?;
  let edge_mutation_counts =
    edge_state_change_counts(&edge_mutations, reconstruction.alphabet()).map_err(OperationError::InferenceFailed)?;
  Ok(Some(TimetreeSequences {
    root_sequence,
    edge_mutations,
    edge_mutation_counts,
  }))
}

pub struct TimetreeSequences {
  pub root_sequence: Seq,
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  pub edge_mutation_counts: BTreeMap<GraphEdgeKey, usize>,
}

pub(crate) fn date_branch_lengths(
  graph: &Graph,
  node_dates: &BTreeMap<GraphNodeKey, Option<f64>>,
) -> BTreeMap<GraphEdgeKey, Option<f64>> {
  graph
    .get_edges()
    .map(|edge| {
      let length = node_dates[&edge.source()]
        .zip(node_dates[&edge.target()])
        .map(|(parent, child)| child - parent);
      (edge.key(), length)
    })
    .collect()
}
