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
use crate::coalescent::population_size::effective_population_size;
use crate::coalescent::skyline::SkylineParams;
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
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
use crate::timetree::coalescent::CoalescentOutput;
use crate::timetree::coalescent_timescale::{
  CoalescentMode, build_coalescent_output, coalescent_mode, coalescent_timescale,
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
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
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
