use crate::branch_lengths::branch_length_or_zero;
use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockFit, ClockFitParams, ClockTree};
use crate::clock::clock_set::ClockSet;
use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
use crate::clock::divergence::root_to_node_divergences;
use crate::clock::reroot::RerootParams;
use crate::clock::rtt::{ClockDateSource, ClockRegressionResult, clock_fit_regression_results};
use crate::coalescent::mean_field::CoalescentPrior;
use crate::coalescent::population_size::effective_population_size;
use crate::error::OperationError;
use crate::optimize::params::BranchLengthMode;
use crate::optimize::partition_rates::fix_equal_rates;
use crate::partition::create::Representation;
use crate::partition::marginal::reconstruction::ReconstructionRef;
use crate::partition::result::{PartitionResult, PartitionSequences, PartitionSource};
use crate::partition::set::build_partition_set;
use crate::progress::{LogSink, StageSink};
use crate::seq::alignment::SequenceSource;
use crate::seq::sink::{SeqSink, SeqTrack};
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
use crate::timetree::inference::runner::timetree_branch_lengths;
use crate::timetree::init::{initialize, timetree_gaussian_clock};
use crate::timetree::optimization::reroot::fit_clock_to_dates;
use crate::timetree::params::{
  FinalSequences, TimeMarginalMode, TimetreeContext, TimetreeParams, compute_effective_time_marginal,
};
use crate::timetree::pre_loop::{PreLoopInputs, PreLoopState, run_pre_loop};
use crate::timetree::refinement_loop::run_refinement_loop;
use crate::timetree::round::{RoundInputs, RoundState, final_marginal_round, infer_final_times};
use crate::{progress_info, progress_warn};
use eyre::{Report, WrapErr};
use itertools::izip;
use log::debug;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::date::DateConstraint;
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
    sources,
    dates,
    branch_lengths,
    names,
  } = input;
  let alignment_length = (!sources.is_empty()).then(|| sources.iter().map(|source| source.length).sum());
  let context = prepare_inputs(params, &graph, alignment_length, dates.as_ref(), log)?;

  cancel.check().map_err(OperationError::classify)?;
  stages.report("Clock regression", 0.1, "");
  let (graph, branch_lengths, clock_fit) =
    estimate_initial_clock(params, &context, graph, branch_lengths, &names, log).map_err(OperationError::classify)?;
  let has_alignment = !sources.is_empty();
  let (sources, branch_model) = initialize_branch_model(params, &graph, &branch_lengths, sources, log)?;
  let pre_loop_inputs = PreLoopInputs {
    params,
    context: &context,
    names: &names,
    has_alignment,
  };
  let pre_loop_state = PreLoopState::new(graph, branch_lengths, branch_model, clock_fit);
  let pre_loop =
    run_pre_loop(&pre_loop_inputs, pre_loop_state, cancel, stages, log).map_err(OperationError::classify)?;

  let initial = initialize(params, &context, &names, pre_loop, log).map_err(OperationError::classify)?;
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
  let final_inputs = FinalInputs {
    params,
    context: &context,
    sources,
  };
  assemble_output(final_inputs, seq_sink, initial.outliers, dates, results, log)
}

struct FinalInputs<'a> {
  params: &'a TimetreeParams,
  context: &'a TimetreeContext,
  sources: Vec<PartitionSource>,
}

fn assemble_output(
  inputs: FinalInputs<'_>,
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
  let partition_results = final_partition_results(
    inputs,
    seq_sink,
    &graph,
    branch_model,
    &branch_lengths,
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
    results: partition_results,
    coalescent: results.coalescent_output,
    branch_lengths,
    date_branch_lengths,
    node_dates,
    bad_branches: time_inference.bad_branches,
    divergences: results.divergences,
    outliers,
    names,
    graph,
  })
}

pub struct TimetreeInput {
  pub graph: Graph,
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  pub sources: Vec<SequenceSource>,
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
  pub results: Vec<PartitionResult>,
  pub coalescent: Option<CoalescentOutput>,
  pub dates: Option<BTreeMap<GraphNodeKey, DateConstraint>>,
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
  alignment_length: Option<usize>,
  dates: Option<&BTreeMap<GraphNodeKey, DateConstraint>>,
  log: &dyn LogSink,
) -> Result<TimetreeContext, OperationError> {
  debug!(
    "Branch length mode: {:?}, Keep root: {}",
    params.branch_length_mode, params.keep_root
  );

  let sequence_length = match alignment_length {
    Some(length) => length,
    None if params.branch_length_mode == BranchLengthMode::Marginal => {
      return Err(OperationError::InvalidInput(make_report!(
        "Alignment required for marginal reconstruction"
      )));
    },
    None => params.sequence_length.filter(|&length| length > 0).ok_or_else(|| {
      OperationError::InvalidParams(make_report!(
        "--sequence-length is required without an alignment: the time inference needs the sequence length"
      ))
    })?,
  };
  let time_marginal = compute_effective_time_marginal(params.time_marginal, params.confidence, log);

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
    sequence_length,
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
  let fit = ClockFitParams {
    outliers: &BTreeSet::new(),
    clock_rate: params.clock_rate,
    keep_root: params.keep_root,
    reroot_params: &reroot_params,
    prev_clock_rate: None,
  };
  let (
    ClockTree {
      graph, branch_lengths, ..
    },
    clock_reroot,
  ) = fit_clock_to_dates(graph, branch_lengths, &context.date_constraints, &fit, names, log)
    .wrap_err("Failed to infer clock model")?;
  Ok((graph, branch_lengths, clock_reroot.into_clock_fit()?))
}

fn initialize_branch_model(
  params: &TimetreeParams,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  sources: Vec<SequenceSource>,
  log: &dyn LogSink,
) -> Result<(Vec<PartitionSource>, BranchModel), OperationError> {
  match params.branch_length_mode {
    BranchLengthMode::Input => {
      progress_info!(log, "Branch length mode: Input - using tree branch lengths");
      Ok((vec![], BranchModel::Input))
    },
    BranchLengthMode::Marginal => {
      progress_info!(
        log,
        "Branch length mode: Marginal - initializing partitions from alignment"
      );
      if sources.is_empty() {
        return Err(OperationError::InvalidInput(make_report!(
          "Alignment required for marginal reconstruction"
        )));
      }
      let (sources, partitions) = build_partition_set(
        Representation::resolve(params.dense),
        params.model,
        graph,
        sources,
        &branch_lengths_or_zero(branch_lengths),
        log,
      )
      .map_err(OperationError::classify)?;
      let partitions =
        fix_equal_rates(graph, partitions, params.partition_rates, branch_lengths).map_err(OperationError::classify)?;
      Ok((sources, BranchModel::Marginal(partitions)))
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

/// Rate uncertainty of the Gaussian timetree at the final rate, or `None` when the sequence length is unknown or the
/// approximation degenerates
fn gaussian_rate_std(inputs: &RoundInputs<'_>, state: &RoundState, log: &dyn LogSink) -> Result<Option<f64>, Report> {
  Ok(timetree_gaussian_clock(inputs, state, log)?.and_then(|clock| clock.rate_std))
}

fn refine_final_times(
  inputs: &RoundInputs<'_>,
  coalescent: &CoalescentSetup,
  timescale: &CoalescentTimescale,
  state: RoundState,
  log: &dyn LogSink,
) -> Result<FinalTimes, Report> {
  let params = inputs.params;
  let estimated_rate_std = if params.clock_rate.is_none() {
    gaussian_rate_std(inputs, &state, log)?
  } else {
    None
  };
  let state = if params.clock_rate.is_none() {
    RoundState {
      clock_model: state
        .clock_model
        .estimated_at_rate(&ClockSet::from_points(&state.clock_points), estimated_rate_std),
      ..state
    }
  } else {
    state
  };
  let rate_std = if params.confidence {
    let rate_std = determine_rate_std(params.clock_std_dev, || Ok(estimated_rate_std))?;
    if rate_std.is_none() {
      progress_warn!(
        log,
        "Date intervals leave out the rate uncertainty: it needs --clock-std-dev, or an estimated rate with the \
         alignment or --sequence-length."
      );
    }
    rate_std
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

fn final_partition_results(
  inputs: FinalInputs<'_>,
  mut seq_sink: Option<&mut dyn SeqSink>,
  graph: &Graph,
  branch_model: BranchModel,
  divergence_branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  final_branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<Vec<PartitionResult>, OperationError> {
  let FinalInputs {
    params,
    context,
    sources,
  } = inputs;
  let BranchModel::Marginal(partitions) = branch_model else {
    if params.include_leaves || params.impute_missing_data {
      progress_warn!(
        log,
        "Ignoring tip-state flags (--include-leaves / --impute-missing-data / --reconstruct-tip-states): \
         no ancestral reconstruction was performed under --branch-length-mode=input"
      );
    }
    return Ok(vec![]);
  };
  let partitions = match context.final_sequences {
    FinalSequences::AfterMarginalUpdate => {
      partitions
        .marginal_update(graph, final_branch_lengths)
        .map_err(OperationError::classify)?
        .0
    },
    FinalSequences::Skip | FinalSequences::FromLastRound => partitions,
  };
  let streams_sequences = !matches!(context.final_sequences, FinalSequences::Skip);
  if streams_sequences && let Some(sink) = seq_sink.as_deref_mut() {
    sink.on_topology(graph).map_err(OperationError::SinkFailed)?;
  }
  let indel_rates = params
    .indel_policy
    .edge_indels(graph, &partitions, divergence_branch_lengths)
    .partition_rates(partitions.len());
  izip!(sources, partitions.iter(), indel_rates)
    .map(|(source, partition, indel_rate)| {
      let sequences = streams_sequences
        .then(|| stream_final_sequences(params, seq_sink.as_deref_mut(), graph, partition, source.number))
        .transpose()?;
      Ok(PartitionResult {
        number: source.number,
        path: source.path,
        alphabet: partition.alphabet().clone(),
        length: source.length,
        mask: source.mask,
        model_name: params.model,
        gtr: Some(partition.gtr().clone()),
        indel_rate,
        log_lh: Some(partition.graph_log_lh(graph).map_err(OperationError::classify)?),
        sequences,
      })
    })
    .collect()
}

fn stream_final_sequences(
  params: &TimetreeParams,
  seq_sink: Option<&mut (dyn SeqSink + '_)>,
  graph: &Graph,
  partition: ReconstructionRef<'_>,
  number: usize,
) -> Result<PartitionSequences, OperationError> {
  let mutations = partition.stream_sequences(
    graph,
    params.impute_missing_data,
    SeqTrack::Nuc(number),
    params.include_leaves,
    params.report_ambiguous,
    seq_sink,
  )?;
  PartitionSequences::new(mutations, partition.alphabet()).map_err(OperationError::InferenceFailed)
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
