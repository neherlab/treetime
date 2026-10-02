use crate::alphabet::alphabet::Alphabet;
use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockFit, ClockTree, ClockVarianceParams};
use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
use crate::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::clock::reroot::RerootParams;
use crate::clock::rtt::{ClockDateSource, ClockRegressionResult, clock_fit_regression_results};
use crate::coalescent::coalescent::CoalescentModel;
use crate::coalescent::population_size::effective_population_size;
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::optimize::params::BranchLengthMode;
use crate::partition::create::{Representation, build_marginal_partition};
use crate::partition::marginal::reconstruction::MarginalReconstruction;
use crate::partition::marginal::sequences::emitted_nodes;
use crate::progress::{LogSink, StageSink};
use crate::seq::alignment::node_seq_inputs;
use crate::seq::sink::{SeqItem, SeqSink, SeqTrack};
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
use crate::timetree::divergence::final_divergences;
use crate::timetree::inference::runner::timetree_branch_lengths;
use crate::timetree::inference::time_inference::TimeInference;
use crate::timetree::optimization::reroot::{DatedClockFit, fit_clock_to_dates};
use crate::timetree::params::{
  TimeMarginalMode, TimetreeContext, TimetreeParams, build_covariation_clock_params, compute_effective_time_marginal,
};
use crate::timetree::pre_loop::{PreLoopInputs, PreLoopState, run_pre_loop};
use crate::timetree::refinement_loop::run_refinement_loop;
use crate::timetree::round::{RoundInputs, RoundState, final_marginal_round, infer_final_times, run_initial_round};
use crate::{progress_info, progress_warn};
use eyre::{Report, WrapErr};
use log::debug;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::AlignmentRecord;
use treetime_primitives::date::DatesMap;
use treetime_utils::make_report;

pub fn run(
  params: &TimetreeParams,
  input: TimetreeInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  trace_sink: Option<Box<dyn TraceSink + '_>>,
  seq_sink: Option<Box<dyn SeqSink>>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<TimetreeOutput, OperationError> {
  progress_info!(log, "# TreeTime Timetree Estimation");
  validate_params(params)?;
  let context = prepare_inputs(params, &input, names, log)?;

  cancel.check()?;
  stages.report("Clock regression", 0.1, "");
  let (graph, branch_lengths, clock_fit) =
    estimate_initial_clock(params, &context, input.graph, input.branch_lengths, names, log)?;
  let aln = input.sequences.as_deref();
  let init = initialize_branch_model(params, &graph, &branch_lengths, input.alphabet, aln, names, log)?;
  let pre_loop_inputs = PreLoopInputs {
    params,
    context: &context,
    names,
    has_alignment: aln.is_some(),
  };
  let pre_loop_state = PreLoopState::new(graph, branch_lengths, init.branch_model, clock_fit);
  let pre_loop = run_pre_loop(&pre_loop_inputs, pre_loop_state, cancel, stages, log)?;

  let initial = run_initial_round(params, &context, names, pre_loop, log)?;
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
  )?;
  report_coalescent_size(params, coalescent.mode, &timescale, log);

  cancel.check()?;
  stages.report("Postprocessing", 0.85, "");
  progress_info!(log, "### TreeTime: postprocessing");
  let final_times = refine_final_times(&round_inputs, coalescent, &timescale, state, log)?;
  let filter_divergences = initial.filter_divergences.as_ref();
  let results = gather_results(
    params,
    &context,
    coalescent,
    &timescale,
    filter_divergences,
    final_times,
  )?;
  let state = emit_sequences(params, context.time_marginal, seq_sink, results.state, log)?;

  Ok(TimetreeOutput {
    graph: state.graph,
    clock_model: state.clock_model,
    clock_regression: results.clock_regression,
    confidence_intervals: results.confidence_intervals,
    partitions: match state.branch_model {
      BranchModel::Input => vec![],
      BranchModel::Marginal(partition) => vec![partition],
    },
    dates: input.dates,
    gtr: init.gtr,
    model_name: init.model_name,
    coalescent: results.coalescent_output,
    rate_susceptibility_dates: results.rate_susceptibility_dates,
    clock_branch_lengths: state.clock_branch_lengths,
    branch_lengths: state.branch_lengths,
    divergences: results.divergences,
    outliers: initial.outliers,
    time_inference: state.time_inference,
    gammas: state.gammas,
    names: state.names,
  })
}

pub struct TimetreeInput {
  pub graph: Graph,
  pub alphabet: Alphabet,
  pub sequences: Option<Vec<AlignmentRecord>>,
  pub dates: Option<DatesMap>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

#[derive(Serialize)]
pub struct TimetreeOutput {
  #[serde(skip)]
  pub graph: Graph,
  #[serde(skip)]
  pub clock_model: ClockModel,
  #[serde(skip)]
  pub clock_regression: Vec<ClockRegressionResult>,
  #[serde(skip)]
  pub confidence_intervals: Option<Vec<NodeConfidenceInterval>>,
  #[serde(skip)]
  pub partitions: Vec<MarginalReconstruction>,
  #[serde(skip)]
  pub dates: Option<DatesMap>,
  #[serde(skip)]
  pub gtr: Option<GTR>,
  #[serde(skip)]
  pub model_name: Option<GtrModelName>,
  #[serde(skip)]
  pub coalescent: Option<CoalescentOutput>,
  #[serde(skip)]
  pub rate_susceptibility_dates: BTreeMap<GraphNodeKey, [f64; 3]>,
  #[serde(skip)]
  pub clock_branch_lengths: BTreeMap<GraphEdgeKey, f64>,
  #[serde(skip)]
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  #[serde(skip)]
  pub divergences: BTreeMap<GraphNodeKey, f64>,
  #[serde(skip)]
  pub outliers: BTreeSet<GraphNodeKey>,
  #[serde(skip)]
  pub time_inference: TimeInference,
  #[serde(skip)]
  pub gammas: BTreeMap<GraphEdgeKey, f64>,
  #[serde(skip)]
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
}

fn validate_params(params: &TimetreeParams) -> Result<(), OperationError> {
  if params.n_branches_posterior.is_some() {
    return Err(OperationError::InvalidParams(make_report!(
      "--n-branches-posterior is not yet implemented"
    )));
  }
  Ok(())
}

fn prepare_inputs(
  params: &TimetreeParams,
  input: &TimetreeInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<TimetreeContext, OperationError> {
  debug!(
    "Branch length mode: {:?}, Keep root: {}",
    params.branch_length_mode, params.keep_root
  );

  let time_marginal = compute_effective_time_marginal(
    params.time_marginal,
    params.confidence,
    params.clock_std_dev,
    params.covariation,
    log,
  );

  let covariation_clock_params = build_covariation_clock_params(
    params.covariation,
    params.sequence_length,
    params.tip_slack,
    input.sequences.as_deref(),
    log,
  )
  .map_err(OperationError::InvalidParams)?;

  let date_constraints = if let Some(dates) = &input.dates {
    load_date_constraints(dates, &input.graph, names, log)
      .wrap_err("Failed to load date constraints")
      .map_err(OperationError::InvalidInput)?
  } else {
    DateConstraints::default()
  };

  Ok(TimetreeContext {
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

struct BranchModelInit {
  branch_model: BranchModel,
  gtr: Option<GTR>,
  model_name: Option<GtrModelName>,
}

fn initialize_branch_model(
  params: &TimetreeParams,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  alphabet: Alphabet,
  aln: Option<&[AlignmentRecord]>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<BranchModelInit, Report> {
  match params.branch_length_mode {
    BranchLengthMode::Input => {
      progress_info!(log, "Branch length mode: Input - using tree branch lengths");
      Ok(BranchModelInit {
        branch_model: BranchModel::Input,
        gtr: None,
        model_name: None,
      })
    },
    BranchLengthMode::Marginal => {
      progress_info!(
        log,
        "Branch length mode: Marginal - initializing partitions from alignment"
      );
      let aln_data = aln.ok_or_else(|| make_report!("Alignment required for marginal reconstruction"))?;
      let node_inputs = node_seq_inputs(graph, names, aln_data.to_vec());
      let reconstruction = build_marginal_partition(
        Representation::resolve(params.dense),
        params.model,
        graph,
        0,
        alphabet,
        &node_inputs,
        &branch_lengths_or_zero(branch_lengths),
        log,
      )?;
      Ok(BranchModelInit {
        gtr: Some(reconstruction.gtr().clone()),
        branch_model: BranchModel::Marginal(reconstruction),
        model_name: Some(params.model),
      })
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

struct FinalTimes {
  state: RoundState,
  rate_std: Option<f64>,
  rate_susceptibility_dates: BTreeMap<GraphNodeKey, [f64; 3]>,
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
    determine_rate_std(params.clock_std_dev, params.covariation, &state.clock_model, log)?
  } else {
    None
  };

  let final_model = CoalescentModel::new(&coalescent.lineage_counts, &timescale.distribution)?;
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

struct FinalResults {
  state: RoundState,
  rate_susceptibility_dates: BTreeMap<GraphNodeKey, [f64; 3]>,
  confidence_intervals: Option<Vec<NodeConfidenceInterval>>,
  coalescent_output: Option<CoalescentOutput>,
  divergences: BTreeMap<GraphNodeKey, f64>,
  clock_regression: Vec<ClockRegressionResult>,
}

fn gather_results(
  params: &TimetreeParams,
  context: &TimetreeContext,
  coalescent: &CoalescentSetup,
  timescale: &CoalescentTimescale,
  filter_divergences: Option<&BTreeMap<GraphNodeKey, f64>>,
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
  });

  let coalescent_output = build_coalescent_output(
    coalescent.mode,
    timescale,
    params.gen_per_year,
    &coalescent.skyline_params,
  )?;

  let divergences = final_divergences(&state.graph, &state.branch_lengths, &state.names, filter_divergences)?;
  let names = assign_node_names(state.names, &state.graph)?;

  let clock_regression = clock_fit_regression_results(&state.clock_model, &state.clock_points, &names, |key| {
    if context.date_constraints.date_constraint(key).is_some() {
      ClockDateSource::Input
    } else {
      ClockDateSource::Inferred
    }
  });

  Ok(FinalResults {
    state: RoundState { names, ..state },
    rate_susceptibility_dates,
    confidence_intervals,
    coalescent_output,
    divergences,
    clock_regression,
  })
}

fn emit_sequences(
  params: &TimetreeParams,
  time_marginal: TimeMarginalMode,
  seq_sink: Option<Box<dyn SeqSink>>,
  state: RoundState,
  log: &dyn LogSink,
) -> Result<RoundState, OperationError> {
  let Some(mut sink) = seq_sink else {
    if matches!(state.branch_model, BranchModel::Input) && (params.include_leaves || params.impute_missing_data) {
      progress_warn!(
        log,
        "Ignoring tip-state flags (--include-leaves / --impute-missing-data / --reconstruct-tip-states): \
         no ancestral reconstruction was performed under --branch-length-mode=input"
      );
    }
    return Ok(state);
  };
  let BranchModel::Marginal(partition) = state.branch_model else {
    return Err(OperationError::InvalidParams(make_report!(
      "Reconstructed sequence output requires ancestral reconstruction; \
       incompatible with --branch-length-mode=input"
    )));
  };
  let graph = &state.graph;

  sink.on_topology(graph).map_err(OperationError::SinkFailed)?;
  let partition = if time_marginal.runs_final_round() {
    partition
  } else {
    partition
      .marginal_update(
        graph,
        &timetree_branch_lengths(graph, &state.branch_lengths, &state.clock_branch_lengths),
      )?
      .0
  };
  for key in emitted_nodes(graph, params.include_leaves)? {
    let seq = partition.node_sequence(graph, params.impute_missing_data, key)?;
    sink.emit(SeqItem {
      key,
      track: SeqTrack::Nuc,
      seq: &seq,
    })?;
  }
  Ok(RoundState {
    branch_model: BranchModel::Marginal(partition),
    ..state
  })
}
