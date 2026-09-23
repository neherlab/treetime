use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::marginal::branch_lengths_or_zero;
use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
use crate::ancestral::sample::SampleMode;
use crate::ancestral::tip_states::TipStates;
use crate::cancel::Cancel;
use crate::clock::clock_filter::clock_filter_inplace;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockVarianceParams, estimate_clock_model_with_reroot_policy};
use crate::clock::clock_state::{ClockInputs, ClockState};
use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
use crate::clock::reroot::RerootParams;
use crate::coalescent::coalescent::CoalescentModel;
use crate::coalescent::lineage_counts::compute_lineage_counts;
use crate::coalescent::node_time::CoalescentNodeTimes;
use crate::coalescent::population_size::effective_population_size;
use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::make_error;
use crate::optimize::dispatch::{run_optimize_mixed, run_optimize_mixed_inner};
use crate::optimize::gather::{
  gather_timetree_edge_contributions, gather_timetree_edge_indel_counts, timetree_total_sequence_length,
};
use crate::optimize::iteration::apply_damping;
use crate::optimize::params::{BranchLengthMode, BranchOptMethod};
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::timetree::marginal::{
  ancestral_reconstruction_timetree, initialize_marginal_timetree, marginal_update_timetree,
};
use crate::partition::timetree::partition::PartitionTimetree;
use crate::progress::ProgressSink;
use crate::seq::alignment::node_seq_inputs;
use crate::seq::sink::{SeqItem, SeqSink, SeqTrack};
use crate::timetree::coalescent::{
  CoalescentBand, CoalescentInputs, CoalescentOutput, CoalescentOutputMode, CoalescentSolve,
};
use crate::timetree::confidence::{
  NodeConfidenceInterval, compute_rate_susceptibility, determine_rate_std, extract_confidence_intervals,
};
use crate::timetree::convergence::optimizer::{IterationContext, TimetreeOptimizer, TraceSink};
use crate::timetree::inference::runner::{commit_clock_branch_lengths, run_timetree, timetree_branch_lengths};
use crate::timetree::optimization::clock_filter::{apply_outlier_bad_branches, report_bad_branches};
use crate::timetree::optimization::reroot::reroot_tree;
use crate::timetree::params::{TimeMarginalMode, build_covariation_clock_params, compute_effective_time_marginal};
use crate::timetree::refinement::{Refinement, RefinementOptions, TopologyRefinement};
use crate::timetree::timetree_state::TimetreeState;
use crate::timetree::utils::initialize_node_divergences;
use eyre::{Report, WrapErr};
use log::{debug, info, warn};
use ndarray::{Array1, array};
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_distribution::Distribution;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_primitives::AlignmentRecord;
use treetime_primitives::date::DatesMap;
use treetime_utils::make_report;
use treetime_utils::sync::random::get_random_number_generator;

const TIMETREE_PRE_STEP_DAMPING: f64 = 0.75;

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub fn run(
  params: &TimetreeParams,
  mut input: TimetreeInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  trace_sink: Option<Box<dyn TraceSink>>,
  mut seq_sink: Option<Box<dyn SeqSink>>,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<TimetreeOutput, OperationError> {
  info!("# TreeTime Timetree Estimation");

  debug!(
    "Branch length mode: {:?}, Keep root: {}",
    params.branch_length_mode, params.keep_root
  );

  let time_marginal = compute_effective_time_marginal(
    params.time_marginal,
    params.confidence,
    params.clock_std_dev,
    params.covariation,
  );

  let covariation_clock_params = build_covariation_clock_params(
    params.covariation,
    params.sequence_length,
    params.tip_slack,
    input.sequences.as_deref(),
  )
  .map_err(OperationError::InvalidParams)?;

  let date_constraints = if let Some(dates) = &input.dates {
    load_date_constraints(dates, &input.graph, names)
      .wrap_err("Failed to load date constraints")
      .map_err(OperationError::InvalidInput)?
  } else {
    DateConstraints::default()
  };

  let mut timetree_state = TimetreeState::seed_from_values(&input.graph, &date_constraints);

  let mut clock_state = ClockState::new(&input.graph);

  let mut branch_lengths = std::mem::take(&mut input.branch_lengths);

  initialize_node_divergences(&input.graph, &mut clock_state, &branch_lengths, names)?;

  let branch_params = BranchPointOptimizationParams::default();

  cancel.check()?;
  progress.report("Clock regression", 0.1, "");
  let reroot_params = RerootParams {
    spec: params.reroot_spec.clone(),
    force_positive_rate: !params.allow_negative_rate,
    ..RerootParams::default()
  };
  clock_state.reseed_transitional(&input.graph);
  let mut clock_inputs = ClockInputs::seed_from_times(&input.graph, &timetree_state.likely_times(&date_constraints));
  let (new_clock_state, clock_reroot) = estimate_clock_model_with_reroot_policy(
    &mut input.graph,
    &mut clock_inputs,
    clock_state,
    &ClockVarianceParams::default(),
    params.clock_rate,
    params.keep_root,
    &branch_params,
    &reroot_params,
    &mut branch_lengths,
    None,
    names,
  )
  .wrap_err("Failed to infer clock model")?;
  clock_state = new_clock_state;
  let mut clock_model = clock_reroot.into_clock_model()?;

  let (mut partitions, partition_gtr, partition_model_name): (
    Vec<PartitionTimetree>,
    Option<GTR>,
    Option<GtrModelName>,
  ) = match params.branch_length_mode {
    BranchLengthMode::Input => {
      info!("Branch length mode: Input - using tree branch lengths");
      (vec![], None, None)
    },
    BranchLengthMode::Marginal => {
      info!("Branch length mode: Marginal - initializing partitions from alignment");
      let init = initialize_partitions_from_params(
        params,
        &input.graph,
        input.alphabet.clone(),
        input.sequences.as_deref(),
        &branch_lengths,
        names,
      )?;
      (init.partitions, Some(init.gtr), Some(init.model_name))
    },
  };

  if let Some(aln) = input.sequences.as_deref() {
    if params.branch_length_mode == BranchLengthMode::Marginal && !partitions.is_empty() {
      info!("### ML branch-length optimization (pre-reroot)");
      let node_inputs = node_seq_inputs(&input.graph, names, aln.to_vec());
      (partitions, _) = initialize_marginal_timetree(
        &input.graph,
        &branch_lengths_or_zero(&branch_lengths),
        partitions,
        &node_inputs,
      )?;
      partitions = optimize_branch_lengths_pre_step(&input.graph, partitions, params.no_indels, &mut branch_lengths)
        .wrap_err("ML branch-length optimization (pre-reroot) failed")?;
    }
  }

  if !params.keep_root {
    info!("First reroot (pre-ancestral)");
    (clock_model, partitions) = reroot_tree(
      &mut input.graph,
      &date_constraints,
      &mut clock_state,
      &timetree_state,
      partitions,
      &ClockVarianceParams::default(),
      params.clock_rate,
      &branch_params,
      &params.reroot_spec,
      !params.allow_negative_rate,
      &mut branch_lengths,
      names,
    )
    .wrap_err("Failed to reroot tree (pre-ancestral)")?;
  }

  if params.clock_filter > 0.0 {
    timetree_state.reseed_from_values(&input.graph);

    let given_dates = timetree_state.likely_times(&date_constraints);
    clock_state.reseed_transitional(&input.graph);
    let clock_inputs = ClockInputs::seed_from_times(&input.graph, &given_dates);
    let result = clock_filter_inplace(
      &input.graph,
      &clock_inputs,
      &mut clock_state,
      &clock_model,
      &branch_lengths,
      params.clock_filter,
    )?;
    report_bad_branches(
      &input.graph,
      &clock_state,
      &clock_model,
      result.iqd,
      &given_dates,
      names,
    );
    apply_outlier_bad_branches(&input.graph, &clock_state, &mut timetree_state)?;
  }

  if let Some(aln) = input.sequences.as_deref() {
    match params.branch_length_mode {
      BranchLengthMode::Input => {
        info!("Using input branch lengths for timetree inference");
      },
      BranchLengthMode::Marginal => {
        info!("### ML branch-length optimization (post-reroot)");
        (partitions, _) = marginal_update_timetree(&input.graph, &branch_lengths_or_zero(&branch_lengths), partitions)?;
        partitions = optimize_branch_lengths_pre_step(&input.graph, partitions, params.no_indels, &mut branch_lengths)
          .wrap_err("ML branch-length optimization (post-reroot) failed")?;
      },
    }
  }

  cancel.check()?;
  progress.report("Initial timetree inference", 0.2, "");
  info!("### TreeTime: initial round");

  let default_clock_params = ClockVarianceParams::default();
  let reroot_clock_params = covariation_clock_params.as_ref().unwrap_or(&default_clock_params);

  if !params.keep_root {
    info!("Reroot (post-ancestral)");
    (clock_model, partitions) = reroot_tree(
      &mut input.graph,
      &date_constraints,
      &mut clock_state,
      &timetree_state,
      partitions,
      reroot_clock_params,
      params.clock_rate,
      &branch_params,
      &params.reroot_spec,
      !params.allow_negative_rate,
      &mut branch_lengths,
      names,
    )
    .wrap_err("Failed to reroot tree (post-ancestral)")?;
  }

  let mut clock_branch_lengths: BTreeMap<GraphEdgeKey, f64> = BTreeMap::new();

  let mut names: BTreeMap<GraphNodeKey, Option<String>> = input
    .graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      (key, names.get(&key).cloned().flatten())
    })
    .collect();

  timetree_state = run_timetree(
    &input.graph,
    &date_constraints,
    &partitions,
    &branch_lengths,
    &names,
    &clock_model,
    None,
    params.no_indels,
    timetree_state,
    &mut clock_state,
  )?;

  let skyline_params = SkylineParams {
    n_points: params.skyline_n_points,
    stiffness: params.skyline_stiffness,
    n_std: params.coalescent_confidence,
    ..SkylineParams::default()
  };

  if params.n_branches_posterior.is_some() {
    return Err(OperationError::InvalidParams(make_report!(
      "--n-branches-posterior is not yet implemented"
    )));
  }
  let coalescent = coalescent_mode(params.coalescent, params.coalescent_opt, params.coalescent_skyline);

  let coalescent_node_times = timetree_state.coalescent_node_times();

  let lineage_counts = compute_lineage_counts(&input.graph, &coalescent_node_times)
    .wrap_err("Failed to compute coalescent lineage counts")?;

  let mut coalescent_tc = coalescent_timescale(coalescent, &input.graph, &skyline_params, &coalescent_node_times)?;

  let prior_wanted = coalescent != CoalescentMode::Disabled;

  if prior_wanted {
    let prior = CoalescentModel::new(&lineage_counts, &coalescent_tc.distribution)?;
    timetree_state = run_timetree(
      &input.graph,
      &date_constraints,
      &partitions,
      &branch_lengths,
      &names,
      &clock_model,
      Some(&prior),
      params.no_indels,
      timetree_state,
      &mut clock_state,
    )?;
  }
  commit_clock_branch_lengths(
    &input.graph,
    clock_model.clock_rate(),
    1.0,
    &mut clock_branch_lengths,
    &timetree_state,
  );

  cancel.check()?;
  progress.report("Optimization", 0.3, "");
  info!("### TreeTime: Optimisation rounds");
  let mut optimizer = TimetreeOptimizer::new(params.max_iter, false);
  if let Some(sink) = trace_sink {
    optimizer = optimizer.with_trace_sink(sink);
  }
  let refinement_options = RefinementOptions {
    relax: params.relax.clone(),
    topology: if params.resolve_polytomies {
      TopologyRefinement::Resolve
    } else {
      TopologyRefinement::Disabled
    },
    clock_rate: params.clock_rate,
    no_indels: params.no_indels,
  };
  let max_iter = params.max_iter;

  let seed = params.seed.unwrap_or_else(rand::random);
  if params.resolve_polytomies {
    info!("Polytomy resolution is stochastic; seed {seed} (pass --seed to reproduce this run)");
  }
  let mut rng = get_random_number_generator(Some(seed));

  while let Some(IterationContext { i }) = optimizer.next_iter() {
    cancel.check()?;
    let iter_fraction = 0.3 + 0.5 * (i as f64 / max_iter as f64);
    progress.report(
      "Optimization",
      iter_fraction,
      &format!("iteration {}/{max_iter}", i + 1),
    );

    if coalescent.is_optimized() {
      coalescent_tc = coalescent_timescale(
        coalescent,
        &input.graph,
        &skyline_params,
        &timetree_state.coalescent_node_times(),
      )?;
    }
    let coalescent_model = CoalescentModel::new(&lineage_counts, &coalescent_tc.distribution)?;
    let merger_rate = coalescent_model.branch_merger_rate_schedule(&coalescent_tc.schedule)?;

    let (new_timetree_state, new_partitions, outcome) = Refinement {
      graph: &mut input.graph,
      partitions,
      clock_model: &mut clock_model,
      clock_params: reroot_clock_params,
      branch_params: &branch_params,
      merger_rate: &merger_rate,
      prior: prior_wanted.then_some(&coalescent_model),
      rng: &mut rng,
      options: &refinement_options,
      constraints: &date_constraints,
      state: timetree_state,
      clock_state: &mut clock_state,
      clock_branch_lengths: &mut clock_branch_lengths,
      branch_lengths: &mut branch_lengths,
      names: &mut names,
    }
    .run()
    .wrap_err_with(|| format!("When running round {i}"))?;
    timetree_state = new_timetree_state;
    partitions = new_partitions;

    optimizer
      .record(
        outcome.sequence_changes,
        outcome.topology.resolved_nodes(),
        outcome.time_change,
        &input.graph,
        &partitions,
        &timetree_state,
        prior_wanted.then_some(&coalescent_tc.distribution),
      )
      .wrap_err("Failed to record convergence metrics")
      .wrap_err_with(|| format!("When running round {i}"))?;
  }

  if coalescent.output_mode().is_some() {
    let tc_values = coalescent_tc.schedule.values();
    info!(
      "Coalescent effective population size (gen_per_year={:.4}, {} segment(s)):",
      params.gen_per_year,
      tc_values.len()
    );
    for (i, &tc) in tc_values.iter().enumerate() {
      info!(
        "  segment {i}: Tc = {tc:.6e}  N_e = {:.6e}",
        effective_population_size(tc, params.gen_per_year)
      );
    }
  }

  cancel.check()?;
  progress.report("Postprocessing", 0.85, "");
  info!("### TreeTime: postprocessing");

  let rate_std = if params.confidence {
    determine_rate_std(params.clock_std_dev, params.covariation, &clock_model)?
  } else {
    None
  };

  let final_model = CoalescentModel::new(&lineage_counts, &coalescent_tc.distribution)?;
  let final_prior = prior_wanted.then_some(&final_model);

  let rate_susceptibility_dates = if let Some(rate_std) = rate_std {
    info!("### Rate susceptibility analysis (rate_std={rate_std:.6e})");
    compute_rate_susceptibility(
      &input.graph,
      &date_constraints,
      &partitions,
      &clock_model,
      final_prior,
      rate_std,
      params.no_indels,
      &branch_lengths,
      &mut timetree_state,
      &mut clock_state,
      &names,
    )
    .wrap_err("Rate susceptibility analysis failed")?
  } else {
    BTreeMap::new()
  };

  if time_marginal == TimeMarginalMode::OnlyFinal {
    info!("### Final round: marginal reconstruction for confidence intervals");
    timetree_state = run_timetree(
      &input.graph,
      &date_constraints,
      &partitions,
      &branch_lengths,
      &names,
      &clock_model,
      final_prior,
      params.no_indels,
      timetree_state,
      &mut clock_state,
    )
    .wrap_err("Final timetree inference failed")?;

    commit_clock_branch_lengths(
      &input.graph,
      clock_model.clock_rate(),
      1.0,
      &mut clock_branch_lengths,
      &timetree_state,
    );

    if !partitions.is_empty() {
      (partitions, _) = marginal_update_timetree(
        &input.graph,
        &timetree_branch_lengths(&input.graph, &branch_lengths, &clock_branch_lengths),
        partitions,
      )?;
    }
  }

  let confidence_intervals = (matches!(time_marginal, TimeMarginalMode::OnlyFinal | TimeMarginalMode::Always)
    || rate_std.is_some())
  .then(|| extract_confidence_intervals(&input.graph, &timetree_state, &rate_susceptibility_dates, &names));

  let coalescent_output = build_coalescent_output(coalescent, &coalescent_tc, params.gen_per_year, &skyline_params)?;

  let names = assign_node_names(names, &input.graph)?;

  if seq_sink.is_some() || params.include_leaves || params.impute_missing_data {
    if partitions.is_empty() {
      if seq_sink.is_some() {
        return Err(OperationError::InvalidParams(make_report!(
          "Reconstructed sequence output requires ancestral reconstruction; \
           incompatible with --branch-length-mode=input"
        )));
      }
      warn!(
        "Ignoring tip-state flags (--include-leaves / --impute-missing-data / --reconstruct-tip-states): \
         no ancestral reconstruction was performed under --branch-length-mode=input"
      );
    } else {
      if let Some(sink) = seq_sink.as_mut() {
        sink.on_topology(&input.graph).map_err(OperationError::SinkFailed)?;
      }
      let branch_lengths_final = timetree_branch_lengths(&input.graph, &branch_lengths, &clock_branch_lengths);
      (partitions, _) = marginal_update_timetree(&input.graph, &branch_lengths_final, partitions)?;
      let mut rng = get_random_number_generator(params.seed);
      ancestral_reconstruction_timetree(
        &input.graph,
        TipStates {
          include_leaves: params.include_leaves,
          impute: params.impute_missing_data,
        },
        &mut partitions,
        SampleMode::Argmax,
        &mut rng,
        |key, seq| match seq_sink.as_mut() {
          Some(sink) => sink.emit(SeqItem {
            key,
            track: SeqTrack::Nuc,
            seq,
          }),
          None => Ok(()),
        },
      )?;
    }
  }

  progress.report("Done", 1.0, "");
  Ok(TimetreeOutput {
    graph: input.graph,
    clock_model,
    confidence_intervals,
    partitions,
    dates: input.dates,
    gtr: partition_gtr,
    model_name: partition_model_name,
    coalescent: coalescent_output,
    rate_susceptibility_dates,
    clock_branch_lengths,
    branch_lengths,
    clock_state,
    timetree_state,
    names,
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
  pub clock_state: ClockState,
  #[serde(skip)]
  pub timetree_state: TimetreeState,
  #[serde(skip)]
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
}

fn coalescent_mode(coalescent: Option<f64>, coalescent_opt: bool, coalescent_skyline: bool) -> CoalescentMode {
  if coalescent_skyline {
    CoalescentMode::Skyline
  } else if coalescent_opt {
    CoalescentMode::Constant
  } else if let Some(tc) = coalescent {
    CoalescentMode::Fixed(tc)
  } else {
    CoalescentMode::Disabled
  }
}

fn coalescent_timescale(
  mode: CoalescentMode,
  graph: &Graph,
  skyline_params: &SkylineParams,
  node_times: &CoalescentNodeTimes,
) -> Result<CoalescentTimescale, Report> {
  let mode = match mode {
    CoalescentMode::Disabled => CoalescentMode::Constant,
    mode => mode,
  };
  estimate_coalescent_tc(mode, graph, skyline_params, node_times)
    .wrap_err("Failed to estimate the coalescent timescale")?
    .ok_or_else(|| make_report!("A coalescent Tc is required, but {mode:?} yielded none"))
}

fn estimate_coalescent_tc(
  mode: CoalescentMode,
  graph: &Graph,
  skyline_params: &SkylineParams,
  node_times: &CoalescentNodeTimes,
) -> Result<Option<CoalescentTimescale>, Report> {
  let n_points = match mode {
    CoalescentMode::Disabled => return Ok(None),
    CoalescentMode::Fixed(tc) => return fixed_timescale(tc, graph, node_times).map(Some),
    CoalescentMode::Constant => 1,
    CoalescentMode::Skyline => skyline_params.n_points,
  };
  let result = optimize_skyline(
    graph,
    &SkylineParams {
      n_points,
      ..skyline_params.clone()
    },
    node_times,
  )?;
  Ok(Some(CoalescentTimescale {
    distribution: result.tc_distribution,
    schedule: result.tc_schedule,
    report: Some(CoalescentTcReport {
      segment_boundaries: result.segment_boundaries,
      band: Some(CoalescentReportBand {
        lower: result.tc_lower_bounds,
        upper: result.tc_upper_bounds,
      }),
      log_likelihood: Some(result.log_likelihood.value()),
    }),
  }))
}

fn fixed_timescale(tc: f64, graph: &Graph, node_times: &CoalescentNodeTimes) -> Result<CoalescentTimescale, Report> {
  let lineage_counts =
    compute_lineage_counts(graph, node_times).wrap_err("Failed to compute coalescent lineage counts")?;
  let breakpoints = lineage_counts.breakpoints();
  if breakpoints.is_empty() {
    return make_error!("Cannot report a fixed coalescent Tc: the tree has no node times to span");
  }
  let t_min = breakpoints[0];
  let t_max = breakpoints[breakpoints.len() - 1];
  Ok(CoalescentTimescale {
    report: Some(CoalescentTcReport {
      segment_boundaries: array![t_min, t_max],
      band: None,
      log_likelihood: None,
    }),
    ..CoalescentTimescale::constant(tc)
  })
}

fn build_coalescent_output(
  requested: CoalescentMode,
  timescale: &CoalescentTimescale,
  gen_per_year: f64,
  skyline_params: &SkylineParams,
) -> Result<Option<CoalescentOutput>, Report> {
  let Some(mode) = requested.output_mode() else {
    return Ok(None);
  };
  let report = timescale.report.as_ref().ok_or_else(|| {
    make_report!("A coalescent output ({mode:?}) must carry a per-segment report, but none was produced")
  })?;

  let tc_values = timescale.schedule.values().to_vec();
  let boundaries = report.segment_boundaries.to_vec();

  let (n_points, stiffness) = match mode {
    CoalescentOutputMode::Skyline => (Some(skyline_params.n_points), Some(skyline_params.stiffness)),
    CoalescentOutputMode::Fixed | CoalescentOutputMode::Constant => (None, None),
  };
  let confidence_n_std = report.band.as_ref().map(|_| skyline_params.n_std);

  let (lower, upper) = match &report.band {
    Some(band) => (band.lower.to_vec(), band.upper.to_vec()),
    None => (Vec::new(), Vec::new()),
  };
  let band = report.band.as_ref().map(|_| CoalescentBand {
    lower: &lower,
    upper: &upper,
  });

  let output = CoalescentOutput::new(
    CoalescentInputs {
      mode,
      n_points,
      stiffness,
      confidence_n_std,
      gen_per_year,
    },
    &CoalescentSolve {
      segment_boundaries: &boundaries,
      tc_values: &tc_values,
      band,
      log_likelihood: report.log_likelihood,
    },
  )?;
  Ok(Some(output))
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum CoalescentMode {
  Disabled,
  Fixed(f64),
  Constant,
  Skyline,
}

impl CoalescentMode {
  fn is_optimized(self) -> bool {
    matches!(self, CoalescentMode::Constant | CoalescentMode::Skyline)
  }

  fn output_mode(self) -> Option<CoalescentOutputMode> {
    match self {
      CoalescentMode::Disabled => None,
      CoalescentMode::Fixed(_) => Some(CoalescentOutputMode::Fixed),
      CoalescentMode::Constant => Some(CoalescentOutputMode::Constant),
      CoalescentMode::Skyline => Some(CoalescentOutputMode::Skyline),
    }
  }
}

struct CoalescentTimescale {
  distribution: Distribution,
  schedule: PiecewiseConstantFn,
  report: Option<CoalescentTcReport>,
}

impl CoalescentTimescale {
  fn constant(tc: f64) -> Self {
    Self {
      distribution: Distribution::constant(tc),
      schedule: PiecewiseConstantFn::new(array![], array![tc]),
      report: None,
    }
  }
}

struct CoalescentTcReport {
  segment_boundaries: Array1<f64>,
  band: Option<CoalescentReportBand>,
  log_likelihood: Option<f64>,
}

struct CoalescentReportBand {
  lower: Array1<f64>,
  upper: Array1<f64>,
}

fn initialize_partitions_from_params(
  params: &TimetreeParams,
  graph: &Graph,
  alphabet: Alphabet,
  aln: Option<&[AlignmentRecord]>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<PartitionInitResult, Report> {
  let model_name = params.model;

  let aln_data = aln.ok_or_else(|| make_report!("Alignment required for marginal reconstruction"))?;
  let node_inputs = node_seq_inputs(graph, names, aln_data.to_vec());
  let created = create_marginal_partition(
    graph,
    0,
    alphabet,
    &node_inputs,
    model_name,
    params.dense,
    &branch_lengths_or_zero(branch_lengths),
  )?;

  let gtr = created.gtr.clone();

  let partition = match created.partition {
    MarginalPartition::Sparse(partition, node_states) => {
      PartitionTimetree::Sparse(SparseReconstruction::seeded(partition, created.gtr, node_states))
    },
    MarginalPartition::Dense(partition) => {
      PartitionTimetree::Dense(DenseReconstruction::seeded(partition, created.gtr, BTreeMap::new()))
    },
  };

  Ok(PartitionInitResult {
    partitions: vec![partition],
    gtr,
    model_name: created.model_name,
  })
}

pub struct TimetreeParams {
  pub model: GtrModelName,
  pub alphabet_name: crate::alphabet::alphabet::AlphabetName,
  pub dense: Option<bool>,
  pub gap_fill: crate::seq::gap_fill::GapFill,
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

struct PartitionInitResult {
  partitions: Vec<PartitionTimetree>,
  gtr: GTR,
  model_name: GtrModelName,
}

fn optimize_branch_lengths_pre_step(
  graph: &Graph,
  partitions: Vec<PartitionTimetree>,
  no_indels: bool,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<Vec<PartitionTimetree>, Report> {
  let old_branch_lengths = branch_lengths.clone();

  {
    let total_length = timetree_total_sequence_length(&partitions);
    let contributions = gather_timetree_edge_contributions(graph, &partitions)?;
    let indel_counts = gather_timetree_edge_indel_counts(graph, &partitions);
    if no_indels {
      run_optimize_mixed_inner(
        graph,
        total_length,
        &contributions,
        &indel_counts,
        BranchOptMethod::BrentSqrt,
        0.0,
        true,
        branch_lengths,
      )
      .wrap_err("ML branch-length optimization pre-step failed")?;
    } else {
      run_optimize_mixed(
        graph,
        total_length,
        &contributions,
        &indel_counts,
        BranchOptMethod::BrentSqrt,
        branch_lengths,
      )
      .wrap_err("ML branch-length optimization pre-step failed")?;
    }
  }

  apply_damping(branch_lengths, &old_branch_lengths, TIMETREE_PRE_STEP_DAMPING, 0);
  let (partitions, _) = marginal_update_timetree(graph, &branch_lengths_or_zero(branch_lengths), partitions)?;

  Ok(partitions)
}

#[cfg(test)]
mod tests {
  use super::{
    CoalescentBand, CoalescentInputs, CoalescentMode, CoalescentOutput, CoalescentOutputMode, CoalescentReportBand,
    CoalescentSolve, CoalescentTcReport, CoalescentTimescale, build_coalescent_output, coalescent_mode,
    estimate_coalescent_tc,
  };
  use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
  use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
  use crate::timetree::timetree_state::TimetreeState;
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_distribution::Distribution;
  use treetime_graph::graph::Graph;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::{o, pretty_assert_array_eq};

  const GEN_PER_YEAR: f64 = 50.0;
  const N_STD: f64 = 2.0;

  #[rustfmt::skip]
  #[rstest]
  #[case::disabled(       None,       false, false, CoalescentMode::Disabled)]
  #[case::fixed(          Some(0.25), false, false, CoalescentMode::Fixed(0.25))]
  #[case::opt_default(    None,       true,  false, CoalescentMode::Constant)]
  #[case::opt_value_ignored(Some(0.25), true, false, CoalescentMode::Constant)]
  #[case::skyline_default(None,       false, true,  CoalescentMode::Skyline)]
  #[case::skyline_over_opt(Some(0.25), true, true,  CoalescentMode::Skyline)]
  #[trace]
  fn test_pipeline_coalescent_mode(
    #[case] coalescent: Option<f64>,
    #[case] coalescent_opt: bool,
    #[case] coalescent_skyline: bool,
    #[case] expected: CoalescentMode,
  ) {
    let actual = coalescent_mode(coalescent, coalescent_opt, coalescent_skyline);

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_pipeline_build_coalescent_output_disabled_returns_none() -> Result<(), Report> {
    let timescale = CoalescentTimescale::constant(1.0);
    let params = SkylineParams {
      n_std: N_STD,
      ..SkylineParams::default()
    };

    let actual = build_coalescent_output(CoalescentMode::Disabled, &timescale, GEN_PER_YEAR, &params)?;

    assert_eq!(None, actual);
    Ok(())
  }

  #[test]
  fn test_pipeline_build_coalescent_output_fixed_emits_one_segment_no_band() -> Result<(), Report> {
    let (graph, constraints) = dated_tree()?;
    let params = SkylineParams {
      n_std: N_STD,
      ..SkylineParams::default()
    };
    let node_times = TimetreeState::seed_from_values(&graph, &constraints).coalescent_node_times();
    let timescale = estimate_coalescent_tc(CoalescentMode::Fixed(2.5), &graph, &params, &node_times)?
      .expect("a fixed Tc yields a coalescent timescale");

    let actual = build_coalescent_output(CoalescentMode::Fixed(2.5), &timescale, GEN_PER_YEAR, &params)?
      .expect("a fixed Tc writes a coalescent output");

    let expected = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Fixed,
        n_points: None,
        stiffness: None,
        confidence_n_std: None,
        gen_per_year: GEN_PER_YEAR,
      },
      &CoalescentSolve {
        segment_boundaries: &[2000.0, 2010.0],
        tc_values: &[2.5],
        band: None,
        log_likelihood: None,
      },
    )?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_pipeline_build_coalescent_output_constant_carries_band() -> Result<(), Report> {
    let params = SkylineParams {
      n_std: N_STD,
      ..SkylineParams::default()
    };
    let timescale = CoalescentTimescale {
      distribution: Distribution::constant(3.0),
      schedule: PiecewiseConstantFn::new(array![], array![3.0]),
      report: Some(CoalescentTcReport {
        segment_boundaries: array![2000.0, 2020.0],
        band: Some(CoalescentReportBand {
          lower: array![2.0],
          upper: array![4.0],
        }),
        log_likelihood: Some(-7.5),
      }),
    };
    let actual = build_coalescent_output(CoalescentMode::Constant, &timescale, GEN_PER_YEAR, &params)?;

    let expected = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Constant,
        n_points: None,
        stiffness: None,
        confidence_n_std: Some(N_STD),
        gen_per_year: GEN_PER_YEAR,
      },
      &CoalescentSolve {
        segment_boundaries: &[2000.0, 2020.0],
        tc_values: &[3.0],
        band: Some(CoalescentBand {
          lower: &[2.0],
          upper: &[4.0],
        }),
        log_likelihood: Some(-7.5),
      },
    )?;
    assert_eq!(Some(expected), actual);
    Ok(())
  }

  #[test]
  fn test_pipeline_build_coalescent_output_skyline_multi_segment_band() -> Result<(), Report> {
    let params = SkylineParams {
      n_points: 2,
      stiffness: 3.0,
      n_std: N_STD,
      ..SkylineParams::default()
    };
    let timescale = CoalescentTimescale {
      distribution: Distribution::constant(3.0),
      schedule: PiecewiseConstantFn::new(array![2010.0], array![3.0, 5.0]),
      report: Some(CoalescentTcReport {
        segment_boundaries: array![2000.0, 2010.0, 2020.0],
        band: Some(CoalescentReportBand {
          lower: array![2.0, 4.0],
          upper: array![4.0, 6.0],
        }),
        log_likelihood: Some(-9.0),
      }),
    };
    let actual = build_coalescent_output(CoalescentMode::Skyline, &timescale, GEN_PER_YEAR, &params)?;

    let expected = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Skyline,
        n_points: Some(2),
        stiffness: Some(3.0),
        confidence_n_std: Some(N_STD),
        gen_per_year: GEN_PER_YEAR,
      },
      &CoalescentSolve {
        segment_boundaries: &[2000.0, 2010.0, 2020.0],
        tc_values: &[3.0, 5.0],
        band: Some(CoalescentBand {
          lower: &[2.0, 4.0],
          upper: &[4.0, 6.0],
        }),
        log_likelihood: Some(-9.0),
      },
    )?;
    assert_eq!(Some(expected), actual);
    Ok(())
  }

  fn dated_tree() -> Result<(Graph, DateConstraints), Report> {
    let dates: DatesMap = btreemap! {
      o!("root") => Some(DateConstraint::exact(2000.0)),
      o!("x")    => Some(DateConstraint::exact(2005.0)),
      o!("a")    => Some(DateConstraint::exact(2010.0)),
      o!("b")    => Some(DateConstraint::exact(2010.0)),
      o!("c")    => Some(DateConstraint::exact(2010.0)),
    };
    let nwk_parsed = nwk_read_str("((a:1,b:1)x:1,c:1)root:0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let constraints = load_date_constraints(&dates, &graph, &names)?;
    Ok((graph, constraints))
  }

  #[test]
  fn test_pipeline_estimate_coalescent_tc_report_carries_the_skyline_solve() -> Result<(), Report> {
    let (graph, constraints) = dated_tree()?;
    let params = SkylineParams {
      n_points: 3,
      ..SkylineParams::default()
    };

    let node_times = TimetreeState::seed_from_values(&graph, &constraints).coalescent_node_times();
    let solve = optimize_skyline(&graph, &params, &node_times)?;
    let timescale = estimate_coalescent_tc(CoalescentMode::Skyline, &graph, &params, &node_times)?
      .expect("skyline mode yields a coalescent timescale");
    let report = timescale
      .report
      .expect("an inferred skyline carries a per-segment report");

    pretty_assert_array_eq!(solve.segment_boundaries, report.segment_boundaries);
    pretty_assert_array_eq!(solve.tc_values, timescale.schedule.values().clone());
    assert_eq!(Some(solve.log_likelihood.value()), report.log_likelihood);
    let band = report.band.expect("an inferred skyline carries a confidence band");
    pretty_assert_array_eq!(solve.tc_lower_bounds, band.lower);
    pretty_assert_array_eq!(solve.tc_upper_bounds, band.upper);

    Ok(())
  }
}
