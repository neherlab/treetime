use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::marginal::branch_lengths_or_zero;
use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
use crate::ancestral::reconstruction::ReconstructedSequences;
use crate::ancestral::sample::SampleMode;
use crate::ancestral::tip_states::TipStates;
use crate::cancel::Cancel;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockFit, ClockTree, ClockVarianceParams};
use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
use crate::clock::reroot::RerootParams;
use crate::clock::rtt::{ClockDateSource, ClockRegressionResult, clock_fit_regression_results};
use crate::coalescent::coalescent::CoalescentModel;
use crate::coalescent::lineage_counts::compute_lineage_counts;
use crate::coalescent::population_size::effective_population_size;
use crate::coalescent::skyline::SkylineParams;
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::optimize::params::BranchLengthMode;
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::timetree::partition::PartitionTimetree;
use crate::progress::ProgressSink;
use crate::seq::alignment::node_seq_inputs;
use crate::seq::gap_fill::GapFill;
use crate::seq::sink::{SeqItem, SeqSink, SeqTrack};
use crate::timetree::branch_model::BranchModel;
use crate::timetree::coalescent::CoalescentOutput;
use crate::timetree::coalescent_timescale::{
  CoalescentMode, CoalescentTimescale, build_coalescent_output, coalescent_mode, coalescent_timescale,
};
use crate::timetree::confidence::{
  NodeConfidenceInterval, RateSusceptibility, compute_rate_susceptibility, determine_rate_std,
  extract_confidence_intervals,
};
use crate::timetree::convergence::optimizer::TraceSink;
use crate::timetree::divergence::final_divergences;
use crate::timetree::inference::bad_branches::undated_leaves;
use crate::timetree::inference::runner::{
  CLOCK_BRANCH_LENGTH_UNDAMPED, blended_clock_branch_lengths, run_timetree, timetree_branch_lengths,
};
use crate::timetree::inference::time_inference::{TimeInference, unit_gammas};
use crate::timetree::optimization::reroot::{DatedClockFit, fit_clock_to_dates};
use crate::timetree::optimization_loop::run_optimization_loop;
use crate::timetree::params::{TimeMarginalMode, build_covariation_clock_params, compute_effective_time_marginal};
use crate::timetree::pre_loop::{PreLoopInputs, PreLoopState, run_pre_loop};
use crate::timetree::round::{RoundInputs, RoundState};
use crate::{progress_info, progress_warn};
use eyre::{Report, WrapErr};
use log::debug;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_primitives::AlignmentRecord;
use treetime_primitives::date::DatesMap;
use treetime_utils::make_report;
use treetime_utils::sync::random::get_random_number_generator;

pub fn run(
  params: &TimetreeParams,
  input: TimetreeInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  trace_sink: Option<Box<dyn TraceSink + '_>>,
  seq_sink: Option<Box<dyn SeqSink>>,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<TimetreeOutput, OperationError> {
  progress_info!(progress, "# TreeTime Timetree Estimation");
  let (context, leaf_bad_branches) = prepare_inputs(params, &input, names, progress)?;

  cancel.check()?;
  progress.report("Clock regression", 0.1, "");
  let (graph, branch_lengths, clock_fit) =
    estimate_initial_clock(params, &context, input.graph, input.branch_lengths, names, progress)?;
  let aln = input.sequences.as_deref();
  let init = initialize_branch_model(params, &graph, &branch_lengths, input.alphabet, aln, names, progress)?;
  let pre_loop_inputs = PreLoopInputs::new(params, &context, names, aln.is_some());
  let pre_loop_state = PreLoopState::new(graph, branch_lengths, init.branch_model, clock_fit, leaf_bad_branches);
  let pre_loop = run_pre_loop(&pre_loop_inputs, pre_loop_state, cancel, progress)?;

  let initial = run_initial_round(params, &context, names, pre_loop, progress)?;
  let round_inputs = RoundInputs::new(params, &context, &initial.leaf_bad_branches, &initial.outliers);
  let coalescent = &initial.coalescent;
  let (state, timescale) = run_optimization_loop(
    params,
    &round_inputs,
    coalescent,
    initial.timescale,
    initial.state,
    trace_sink,
    cancel,
    progress,
  )?;
  report_coalescent_size(params, coalescent.mode, &timescale, progress);

  cancel.check()?;
  progress.report("Postprocessing", 0.85, "");
  progress_info!(progress, "### TreeTime: postprocessing");
  let leaf_bad_branches = &initial.leaf_bad_branches;
  let final_times = refine_final_times(
    params,
    &context,
    leaf_bad_branches,
    coalescent,
    &timescale,
    state,
    progress,
  )?;
  let filter_divergences = initial.filter_divergences.as_ref();
  let results = gather_results(
    params,
    &context,
    coalescent,
    &timescale,
    filter_divergences,
    final_times,
  )?;
  let state = emit_sequences(params, seq_sink, results.state, progress)?;

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
  pub partitions: Vec<PartitionTimetree>,
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

pub struct TimetreeParams {
  pub model: GtrModelName,
  pub alphabet_name: AlphabetName,
  pub dense: Option<bool>,
  pub gap_fill: GapFill,
  pub branch_length_mode: BranchLengthMode,
  pub no_indels: bool,
  pub sequence_length: Option<usize>,
  pub clock_rate: Option<f64>,
  pub clock_std_dev: Option<f64>,
  pub keep_root: bool,
  pub reroot_spec: RerootSpec,
  pub allow_negative_rate: bool,
  pub clock_filter: f64,
  pub covariation: bool,
  pub tip_slack: Option<f64>,
  pub max_iter: usize,
  pub resolve_polytomies: bool,
  pub keep_polytomies: bool,
  pub relax: Vec<f64>,
  pub coalescent: Option<f64>,
  pub coalescent_opt: bool,
  pub coalescent_skyline: bool,
  pub skyline_n_points: usize,
  pub skyline_stiffness: f64,
  pub coalescent_confidence: f64,
  pub gen_per_year: f64,
  pub n_branches_posterior: Option<usize>,
  pub time_marginal: TimeMarginalMode,
  pub confidence: bool,
  pub include_leaves: bool,
  pub impute_missing_data: bool,
  pub report_ambiguous: bool,
  pub zero_based: bool,
  pub seed: Option<u64>,
}

pub(crate) struct TimetreeContext {
  pub time_marginal: TimeMarginalMode,
  pub date_constraints: DateConstraints,
  pub covariation_clock_params: ClockVarianceParams,
  pub branch_params: BranchPointOptimizationParams,
}

pub(crate) struct CoalescentSetup {
  pub mode: CoalescentMode,
  pub skyline_params: SkylineParams,
  pub lineage_counts: PiecewiseConstantFn,
}

impl CoalescentSetup {
  pub(crate) fn prior_wanted(&self) -> bool {
    self.mode != CoalescentMode::Disabled
  }
}

fn prepare_inputs(
  params: &TimetreeParams,
  input: &TimetreeInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  progress: &dyn ProgressSink,
) -> Result<(TimetreeContext, BTreeMap<GraphNodeKey, bool>), OperationError> {
  debug!(
    "Branch length mode: {:?}, Keep root: {}",
    params.branch_length_mode, params.keep_root
  );

  let time_marginal = compute_effective_time_marginal(
    params.time_marginal,
    params.confidence,
    params.clock_std_dev,
    params.covariation,
    progress,
  );

  let covariation_clock_params = build_covariation_clock_params(
    params.covariation,
    params.sequence_length,
    params.tip_slack,
    input.sequences.as_deref(),
    progress,
  )
  .map_err(OperationError::InvalidParams)?;

  let date_constraints = if let Some(dates) = &input.dates {
    load_date_constraints(dates, &input.graph, names, progress)
      .wrap_err("Failed to load date constraints")
      .map_err(OperationError::InvalidInput)?
  } else {
    DateConstraints::default()
  };

  let leaf_bad_branches = undated_leaves(&input.graph, &date_constraints);

  let context = TimetreeContext {
    time_marginal,
    date_constraints,
    covariation_clock_params: covariation_clock_params.unwrap_or_default(),
    branch_params: BranchPointOptimizationParams::default(),
  };
  Ok((context, leaf_bad_branches))
}

fn estimate_initial_clock(
  params: &TimetreeParams,
  context: &TimetreeContext,
  graph: Graph,
  branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  progress: &dyn ProgressSink,
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
  ) = fit_clock_to_dates(graph, branch_lengths, &context.date_constraints, &fit, names, progress)?;
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
  progress: &dyn ProgressSink,
) -> Result<BranchModelInit, Report> {
  match params.branch_length_mode {
    BranchLengthMode::Input => {
      progress_info!(progress, "Branch length mode: Input - using tree branch lengths");
      Ok(BranchModelInit {
        branch_model: BranchModel::Input,
        gtr: None,
        model_name: None,
      })
    },
    BranchLengthMode::Marginal => {
      progress_info!(
        progress,
        "Branch length mode: Marginal - initializing partitions from alignment"
      );
      let aln_data = aln.ok_or_else(|| make_report!("Alignment required for marginal reconstruction"))?;
      let node_inputs = node_seq_inputs(graph, names, aln_data.to_vec());
      let created = create_marginal_partition(
        graph,
        0,
        alphabet,
        &node_inputs,
        params.model,
        params.dense,
        &branch_lengths_or_zero(branch_lengths),
        progress,
      )?;
      let gtr = created.gtr.clone();
      let partition = match created.partition {
        MarginalPartition::Sparse(partition, node_states) => {
          PartitionTimetree::Sparse(SparseReconstruction::seeded(partition, created.gtr, node_states))
        },
        MarginalPartition::Dense(partition) => {
          PartitionTimetree::Dense(DenseReconstruction::seeded(partition, created.gtr))
        },
      };
      Ok(BranchModelInit {
        branch_model: BranchModel::Marginal(partition),
        gtr: Some(gtr),
        model_name: Some(created.model_name),
      })
    },
  }
}

struct InitialRound {
  state: RoundState,
  leaf_bad_branches: BTreeMap<GraphNodeKey, bool>,
  outliers: BTreeSet<GraphNodeKey>,
  filter_divergences: Option<BTreeMap<GraphNodeKey, f64>>,
  coalescent: CoalescentSetup,
  timescale: CoalescentTimescale,
}

fn run_initial_round(
  params: &TimetreeParams,
  context: &TimetreeContext,
  input_names: &BTreeMap<GraphNodeKey, Option<String>>,
  pre_loop: PreLoopState,
  progress: &dyn ProgressSink,
) -> Result<InitialRound, OperationError> {
  let PreLoopState {
    graph,
    branch_lengths,
    branch_model,
    clock_fit: ClockFit {
      model: clock_model,
      points: clock_points,
    },
    leaf_bad_branches,
    outliers,
    filter_divergences,
  } = pre_loop;

  let names: BTreeMap<GraphNodeKey, Option<String>> = graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      (key, input_names.get(&key).cloned().flatten())
    })
    .collect();
  let gammas = unit_gammas(&graph);

  let run = |prior: Option<&CoalescentModel>| {
    run_timetree(
      &graph,
      &context.date_constraints,
      &leaf_bad_branches,
      &gammas,
      &branch_model,
      &branch_lengths,
      &names,
      &clock_model,
      prior,
      params.no_indels,
      progress,
    )
  };
  let time_inference = run(None)?;

  if params.n_branches_posterior.is_some() {
    return Err(OperationError::InvalidParams(make_report!(
      "--n-branches-posterior is not yet implemented"
    )));
  }

  let (coalescent, timescale) = setup_coalescent(params, &graph, &time_inference, &names, progress)?;
  let time_inference = if coalescent.prior_wanted() {
    let prior = CoalescentModel::new(&coalescent.lineage_counts, &timescale.distribution)?;
    run(Some(&prior))?
  } else {
    time_inference
  };
  let clock_branch_lengths = blended_clock_branch_lengths(
    &graph,
    clock_model.clock_rate(),
    CLOCK_BRANCH_LENGTH_UNDAMPED,
    &BTreeMap::new(),
    &time_inference.node_times(),
    &gammas,
    progress,
  );

  Ok(InitialRound {
    state: RoundState {
      graph,
      names,
      branch_model,
      branch_lengths,
      clock_model,
      clock_points,
      clock_branch_lengths,
      gammas,
      time_inference,
    },
    leaf_bad_branches,
    outliers,
    filter_divergences,
    coalescent,
    timescale,
  })
}

fn setup_coalescent(
  params: &TimetreeParams,
  graph: &Graph,
  time_inference: &TimeInference,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  progress: &dyn ProgressSink,
) -> Result<(CoalescentSetup, CoalescentTimescale), Report> {
  let skyline_params = SkylineParams {
    n_points: params.skyline_n_points,
    stiffness: params.skyline_stiffness,
    n_std: params.coalescent_confidence,
    ..SkylineParams::default()
  };
  let mode = coalescent_mode(params.coalescent, params.coalescent_opt, params.coalescent_skyline);
  let coalescent_node_times = time_inference.coalescent_node_times()?;
  let lineage_counts =
    compute_lineage_counts(graph, &coalescent_node_times).wrap_err("Failed to compute coalescent lineage counts")?;
  let timescale = coalescent_timescale(mode, graph, &skyline_params, &coalescent_node_times, names, progress)?;
  let setup = CoalescentSetup {
    mode,
    skyline_params,
    lineage_counts,
  };
  Ok((setup, timescale))
}

fn report_coalescent_size(
  params: &TimetreeParams,
  mode: CoalescentMode,
  timescale: &CoalescentTimescale,
  progress: &dyn ProgressSink,
) {
  if mode.output_mode().is_none() {
    return;
  }
  let tc_values = timescale.schedule.values();
  progress_info!(
    progress,
    "Coalescent effective population size (gen_per_year={:.4}, {} segment(s)):",
    params.gen_per_year,
    tc_values.len()
  );
  for (i, &tc) in tc_values.iter().enumerate() {
    progress_info!(
      progress,
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
  params: &TimetreeParams,
  context: &TimetreeContext,
  leaf_bad_branches: &BTreeMap<GraphNodeKey, bool>,
  coalescent: &CoalescentSetup,
  timescale: &CoalescentTimescale,
  state: RoundState,
  progress: &dyn ProgressSink,
) -> Result<FinalTimes, Report> {
  let constraints = &context.date_constraints;
  let rate_std = if params.confidence {
    determine_rate_std(params.clock_std_dev, params.covariation, &state.clock_model, progress)?
  } else {
    None
  };

  let final_model = CoalescentModel::new(&coalescent.lineage_counts, &timescale.distribution)?;
  let final_prior = coalescent.prior_wanted().then_some(&final_model);

  let (state, rate_susceptibility_dates) = if let Some(rate_std) = rate_std {
    progress_info!(progress, "### Rate susceptibility analysis (rate_std={rate_std:.6e})");
    let RateSusceptibility { dates, central } = compute_rate_susceptibility(
      &state.graph,
      constraints,
      leaf_bad_branches,
      &state.gammas,
      &state.branch_model,
      &state.clock_model,
      final_prior,
      rate_std,
      params.no_indels,
      &state.branch_lengths,
      &state.names,
      progress,
    )
    .wrap_err("Rate susceptibility analysis failed")?;
    let state = RoundState {
      time_inference: central,
      ..state
    };
    (state, dates)
  } else {
    (state, BTreeMap::new())
  };

  let state = if context.time_marginal == TimeMarginalMode::OnlyFinal {
    progress_info!(
      progress,
      "### Final round: marginal reconstruction for confidence intervals"
    );
    final_marginal_round(params, constraints, leaf_bad_branches, final_prior, state, progress)?
  } else {
    state
  };

  Ok(FinalTimes {
    state,
    rate_std,
    rate_susceptibility_dates,
  })
}

fn final_marginal_round(
  params: &TimetreeParams,
  constraints: &DateConstraints,
  leaf_bad_branches: &BTreeMap<GraphNodeKey, bool>,
  prior: Option<&CoalescentModel>,
  state: RoundState,
  progress: &dyn ProgressSink,
) -> Result<RoundState, Report> {
  let time_inference = run_timetree(
    &state.graph,
    constraints,
    leaf_bad_branches,
    &state.gammas,
    &state.branch_model,
    &state.branch_lengths,
    &state.names,
    &state.clock_model,
    prior,
    params.no_indels,
    progress,
  )
  .wrap_err("Final timetree inference failed")?;

  let clock_branch_lengths = blended_clock_branch_lengths(
    &state.graph,
    state.clock_model.clock_rate(),
    CLOCK_BRANCH_LENGTH_UNDAMPED,
    &state.clock_branch_lengths,
    &time_inference.node_times(),
    &state.gammas,
    progress,
  );

  let timetree_lengths = timetree_branch_lengths(&state.graph, &state.branch_lengths, &clock_branch_lengths);
  let branch_model = state.branch_model.marginal_update(&state.graph, &timetree_lengths)?;

  Ok(RoundState {
    branch_model,
    clock_branch_lengths,
    time_inference,
    ..state
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
  mut seq_sink: Option<Box<dyn SeqSink>>,
  state: RoundState,
  progress: &dyn ProgressSink,
) -> Result<RoundState, OperationError> {
  if seq_sink.is_none() && !params.include_leaves && !params.impute_missing_data {
    return Ok(state);
  }
  let partition = match state.branch_model {
    BranchModel::Marginal(partition) => partition,
    BranchModel::Input => {
      if seq_sink.is_some() {
        return Err(OperationError::InvalidParams(make_report!(
          "Reconstructed sequence output requires ancestral reconstruction; \
           incompatible with --branch-length-mode=input"
        )));
      }
      progress_warn!(
        progress,
        "Ignoring tip-state flags (--include-leaves / --impute-missing-data / --reconstruct-tip-states): \
         no ancestral reconstruction was performed under --branch-length-mode=input"
      );
      return Ok(RoundState {
        branch_model: BranchModel::Input,
        ..state
      });
    },
  };
  let graph = &state.graph;

  if let Some(sink) = seq_sink.as_mut() {
    sink.on_topology(graph).map_err(OperationError::SinkFailed)?;
  }
  let partition = partition.marginal_update(
    graph,
    &timetree_branch_lengths(graph, &state.branch_lengths, &state.clock_branch_lengths),
  )?;
  let mut rng = get_random_number_generator(params.seed);
  let ReconstructedSequences {
    sequences,
    emitted_nodes,
  } = partition.reconstruct_sequences(
    graph,
    TipStates {
      include_leaves: params.include_leaves,
      impute: params.impute_missing_data,
    },
    SampleMode::Argmax,
    &mut rng,
  )?;
  if let Some(sink) = seq_sink.as_mut() {
    for key in emitted_nodes {
      sink.emit(SeqItem {
        key,
        track: SeqTrack::Nuc,
        seq: &sequences[&key],
      })?;
    }
  }
  Ok(RoundState {
    branch_model: BranchModel::Marginal(partition),
    ..state
  })
}
