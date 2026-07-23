use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::clock::clock_filter::clock_filter_inplace;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
use crate::clock::date_constraints::load_date_constraints;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
use crate::clock::reroot::RerootParams;
use crate::coalescent::optimize_tc::optimize_tc;
use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::make_error;
use crate::optimize::dispatch::{run_optimize_mixed, run_optimize_mixed_inner};
use crate::optimize::iteration::{apply_damping, save_branch_lengths};
use crate::optimize::params::{BranchLengthMode, BranchOptMethod};
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::timetree::{GraphTimetree, PartitionTimetree, PartitionTimetreeAllVec, PartitionTimetreeRef};
use crate::partition::traits::HasGtr;
use crate::progress::ProgressSink;
use crate::timetree::confidence::{
  NodeConfidenceInterval, compute_rate_susceptibility, determine_rate_std, extract_confidence_intervals,
};
use crate::timetree::convergence::optimizer::{IterationContext, TimetreeOptimizer};
use crate::timetree::inference::runner::run_timetree;
use crate::timetree::optimization::clock_filter::{apply_outlier_bad_branches, report_bad_branches};
use crate::timetree::optimization::reroot::reroot_tree;
use crate::timetree::params::{TimeMarginalMode, build_covariation_clock_params, compute_effective_time_marginal};
use crate::timetree::refinement::{RefinementParams, run_refinement_iteration};
use crate::timetree::utils::{initialize_clock_totals_from_time_distributions, initialize_node_divergences};
use eyre::{Report, WrapErr};
use log::{debug, info, warn};
use parking_lot::RwLock;
use serde::Serialize;
use std::io::Write;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_io::dates_csv::DatesMap;
use treetime_io::fasta::FastaRecord;
use treetime_utils::make_report;

const TIMETREE_PRE_STEP_DAMPING: f64 = 0.75;

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
  pub n_skyline: usize,
  /// Skyline smoothing stiffness (penalizes changes in 1/Tc; units of time^2).
  pub skyline_stiffness: f64,
  pub n_branches_posterior: Option<usize>,
  pub time_marginal: TimeMarginalMode,
  pub confidence: bool,
  pub reconstruct_tip_states: bool,
  pub report_ambiguous: bool,
  pub zero_based: bool,
}

pub struct TimetreeInput {
  pub graph: GraphTimetree,
  pub alphabet: Alphabet,
  pub sequences: Option<Vec<FastaRecord>>,
  pub dates: Option<DatesMap>,
}

#[derive(Serialize)]
pub struct TimetreeOutput {
  #[serde(skip)]
  pub graph: GraphTimetree,
  #[serde(skip)]
  pub clock_model: ClockModel,
  #[serde(skip)]
  pub confidence_intervals: Option<Vec<NodeConfidenceInterval>>,
  #[serde(skip)]
  pub partitions: PartitionTimetreeAllVec,
  #[serde(skip)]
  pub dates: Option<DatesMap>,
  #[serde(skip)]
  pub gtr: Option<GTR>,
  #[serde(skip)]
  pub model_name: Option<GtrModelName>,
}

pub fn run(
  params: &TimetreeParams,
  mut input: TimetreeInput,
  tracelog: Option<Box<dyn Write + Send>>,
  progress: &dyn ProgressSink,
) -> Result<TimetreeOutput, Report> {
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
  )?;

  if let Some(dates) = &input.dates {
    load_date_constraints(dates, &input.graph).wrap_err("Failed to load date constraints")?;
  }

  initialize_node_divergences(&input.graph)?;

  let branch_params = BranchPointOptimizationParams::default();

  progress.check_cancelled()?;
  progress.report("Clock regression", 0.1, "");
  let reroot_params = RerootParams {
    spec: params.reroot_spec.clone(),
    force_positive_rate: !params.allow_negative_rate,
    ..RerootParams::default()
  };
  let mut clock_model = estimate_clock_model_with_reroot_policy(
    &mut input.graph,
    &ClockParams::default(),
    params.clock_rate,
    params.keep_root,
    &branch_params,
    &reroot_params,
    None,
  )
  .wrap_err("Failed to infer clock model")?
  .into_clock_model()?;

  let (partitions, partition_gtr, partition_model_name): (PartitionTimetreeAllVec, Option<GTR>, Option<GtrModelName>) =
    match params.branch_length_mode {
      BranchLengthMode::Input => {
        info!("Branch length mode: Input - using tree branch lengths");
        (vec![], None, None)
      },
      BranchLengthMode::Marginal => {
        info!("Branch length mode: Marginal - initializing partitions from alignment");
        let init =
          initialize_partitions_from_params(params, &input.graph, input.alphabet.clone(), input.sequences.as_deref())?;
        (init.partitions, Some(init.gtr), Some(init.model_name))
      },
    };

  if let Some(aln) = input.sequences.as_deref() {
    if params.branch_length_mode == BranchLengthMode::Marginal && !partitions.is_empty() {
      info!("### ML branch-length optimization (pre-reroot)");
      initialize_marginal(&input.graph, &partitions, aln)?;
      optimize_branch_lengths_pre_step(&input.graph, &partitions, params.no_indels)
        .wrap_err("ML branch-length optimization (pre-reroot) failed")?;
    }
  }

  if !params.keep_root {
    info!("First reroot (pre-ancestral)");
    clock_model = reroot_tree(
      &mut input.graph,
      &partitions,
      &ClockParams::default(),
      params.clock_rate,
      &branch_params,
      &params.reroot_spec,
      !params.allow_negative_rate,
    )
    .wrap_err("Failed to reroot tree (pre-ancestral)")?;
  }

  if params.clock_filter > 0.0 {
    let result = clock_filter_inplace(&input.graph, &clock_model, params.clock_filter)?;
    report_bad_branches(&input.graph, &clock_model, result.iqd);
    apply_outlier_bad_branches(&input.graph)?;
  }

  if let Some(aln) = input.sequences.as_deref() {
    match params.branch_length_mode {
      BranchLengthMode::Input => {
        info!("Using input branch lengths for timetree inference");
      },
      BranchLengthMode::Marginal => {
        info!("### ML branch-length optimization (post-reroot)");
        update_marginal(&input.graph, &partitions)?;
        optimize_branch_lengths_pre_step(&input.graph, &partitions, params.no_indels)
          .wrap_err("ML branch-length optimization (post-reroot) failed")?;
      },
    }
  }

  progress.check_cancelled()?;
  progress.report("Initial timetree inference", 0.2, "");
  info!("### TreeTime: initial round");
  info!("### Initializing node times from date constraints");
  initialize_clock_totals_from_time_distributions(&input.graph)?;

  let skyline_params = SkylineParams {
    n_points: params.n_skyline,
    stiffness: params.skyline_stiffness,
    ..SkylineParams::default()
  };

  if params.n_branches_posterior.is_some() {
    return make_error!("--n-branches-posterior is not yet implemented");
  }

  run_timetree(&mut input.graph, &partitions, &clock_model, None, params.no_indels)?;

  let coalescent = coalescent_mode(params.coalescent, params.coalescent_opt, params.coalescent_skyline);

  // Estimate the desired coalescent Tc (constant or skyline) from the first
  // timetree — the earliest point where node times give lineage counts — then
  // re-infer node times under that prior. There is no need to start from a
  // constant Tc and switch to the skyline later: both are cheap analytic solves.
  let mut coalescent_tc = estimate_coalescent_tc(&coalescent, &input.graph, &skyline_params, None)?;

  if coalescent_tc.is_some() {
    run_timetree(
      &mut input.graph,
      &partitions,
      &clock_model,
      coalescent_tc.as_ref(),
      params.no_indels,
    )?;
  }

  let default_clock_params = ClockParams::default();
  let reroot_clock_params = covariation_clock_params.as_ref().unwrap_or(&default_clock_params);

  if !params.keep_root {
    info!("Reroot (post-ancestral)");
    clock_model = reroot_tree(
      &mut input.graph,
      &partitions,
      reroot_clock_params,
      params.clock_rate,
      &branch_params,
      &params.reroot_spec,
      !params.allow_negative_rate,
    )
    .wrap_err("Failed to reroot tree (post-ancestral)")?;

    // Re-establish node times on the rerooted topology before the optimization loop.
    // The reroot can leave internal-node time distributions cleared for the new
    // topology, and the loop's first refinement builds the coalescent model (which
    // reads node times) before its own backward/forward passes recompute them. Run
    // without the coalescent prior: this pass only restores node times, mirroring v0,
    // which recomputes the time tree after rerooting and only adds the merger model
    // inside the iteration loop.
    // packages/legacy/treetime/treetime/treetime.py#L279-L285
    if coalescent_tc.is_some() {
      run_timetree(&mut input.graph, &partitions, &clock_model, None, params.no_indels)
        .wrap_err("Post-reroot timetree inference failed")?;
    }
  }

  progress.check_cancelled()?;
  progress.report("Optimization", 0.3, "");
  info!("### TreeTime: Optimisation rounds");
  let mut optimizer = TimetreeOptimizer::new(params.max_iter, false);
  if let Some(writer) = tracelog {
    optimizer = optimizer.with_tracelog(writer)?;
  }
  let max_iter = params.max_iter;
  while let Some(IterationContext { i }) = optimizer.next_iter() {
    progress.check_cancelled()?;
    let iter_fraction = 0.3 + 0.5 * (i as f64 / max_iter as f64);
    progress.report(
      "Optimization",
      iter_fraction,
      &format!("iteration {}/{max_iter}", i + 1),
    );

    if coalescent.is_optimized() && i >= 2 {
      coalescent_tc = estimate_coalescent_tc(&coalescent, &input.graph, &skyline_params, coalescent_tc.as_ref())?;
    }

    let refinement_params = RefinementParams {
      relax: params.relax.clone(),
      resolve_polytomies: params.resolve_polytomies,
      clock_rate: params.clock_rate,
      no_indels: params.no_indels,
    };
    let (n_diff, n_resolved) = run_refinement_iteration(
      &refinement_params,
      &mut input.graph,
      &partitions,
      &mut clock_model,
      reroot_clock_params,
      &branch_params,
      coalescent_tc.as_ref(),
    )
    .wrap_err_with(|| format!("When running round {i}"))?;

    optimizer
      .record(n_diff, n_resolved, &input.graph, &partitions, coalescent_tc.as_ref())
      .wrap_err("Failed to record convergence metrics")
      .wrap_err_with(|| format!("When running round {i}"))?;
  }

  progress.check_cancelled()?;
  progress.report("Postprocessing", 0.85, "");
  info!("### TreeTime: postprocessing");

  let rate_std = if params.confidence {
    determine_rate_std(params.clock_std_dev, params.covariation, &clock_model)?
  } else {
    None
  };

  if let Some(rate_std) = rate_std {
    info!("### Rate susceptibility analysis (rate_std={rate_std:.6e})");
    compute_rate_susceptibility(
      &mut input.graph,
      &partitions,
      &clock_model,
      coalescent_tc.as_ref(),
      rate_std,
      params.no_indels,
    )
    .wrap_err("Rate susceptibility analysis failed")?;
  }

  if time_marginal == TimeMarginalMode::OnlyFinal {
    info!("### Final round: marginal reconstruction for confidence intervals");
    run_timetree(
      &mut input.graph,
      &partitions,
      &clock_model,
      coalescent_tc.as_ref(),
      params.no_indels,
    )
    .wrap_err("Final timetree inference failed")?;

    if !partitions.is_empty() {
      update_marginal(&input.graph, &partitions)?;
    }
  }

  let confidence_intervals = (matches!(time_marginal, TimeMarginalMode::OnlyFinal | TimeMarginalMode::Always)
    || rate_std.is_some())
  .then(|| extract_confidence_intervals(&input.graph));

  progress.report("Done", 1.0, "");
  Ok(TimetreeOutput {
    graph: input.graph,
    clock_model,
    confidence_intervals,
    partitions,
    dates: input.dates,
    gtr: partition_gtr,
    model_name: partition_model_name,
  })
}

/// The coalescent Tc behavior requested for a run.
#[derive(Clone, Copy, Debug, PartialEq)]
enum CoalescentMode {
  /// No coalescent prior.
  Disabled,
  /// Fixed, user-supplied Tc.
  Fixed(f64),
  /// Optimize a constant Tc each round. `seed` is an optional user-supplied Tc
  /// used only as the fallback if the (analytic) optimization ever fails.
  Constant { seed: Option<f64> },
  /// Optimize a piecewise-constant skyline Tc(t) each round.
  Skyline,
}

impl CoalescentMode {
  /// Whether Tc is re-estimated from the tree, as opposed to fixed or disabled.
  fn is_optimized(self) -> bool {
    matches!(self, CoalescentMode::Constant { .. } | CoalescentMode::Skyline)
  }
}

fn coalescent_mode(coalescent: Option<f64>, coalescent_opt: bool, coalescent_skyline: bool) -> CoalescentMode {
  if coalescent_skyline {
    CoalescentMode::Skyline
  } else if coalescent_opt {
    CoalescentMode::Constant { seed: coalescent }
  } else if let Some(tc) = coalescent {
    CoalescentMode::Fixed(tc)
  } else {
    CoalescentMode::Disabled
  }
}

/// Estimates the coalescent Tc prior for the current tree under `mode`.
///
/// `current` is the previous round's Tc, used as the constant-optimization
/// seed/fallback and preserved when an optimization fails. Skyline and constant
/// optimizations are both cheap analytic solves, so this runs every round rather
/// than deferring the skyline to a final pass.
fn estimate_coalescent_tc(
  mode: &CoalescentMode,
  graph: &GraphTimetree,
  skyline_params: &SkylineParams,
  current: Option<&Distribution>,
) -> Result<Option<Distribution>, Report> {
  match mode {
    CoalescentMode::Disabled => Ok(None),
    CoalescentMode::Fixed(tc) => Ok(Some(Distribution::constant(*tc))),
    CoalescentMode::Constant { seed } => match optimize_tc(graph) {
      Ok(result) => {
        info!(
          "Coalescent Tc = {:.6e} (likelihood = {:.4})",
          result.tc, result.likelihood
        );
        Ok(Some(Distribution::constant(result.tc)))
      },
      // On a degenerate tree prefer the previous round's Tc, then the user seed,
      // then no prior at all (rather than an invented timescale).
      Err(e) => {
        let fallback = current.cloned().or_else(|| (*seed).map(Distribution::constant));
        report_coalescent_failure("Coalescent Tc optimization", &e, fallback.is_some());
        Ok(fallback)
      },
    },
    CoalescentMode::Skyline => match optimize_skyline(graph, skyline_params) {
      Ok(result) => Ok(Some(result.tc_distribution)),
      Err(e) => {
        report_coalescent_failure("Skyline coalescent optimization", &e, current.is_some());
        Ok(current.cloned())
      },
    },
  }
}

/// Reports a coalescent optimization failure and states the resulting prior.
///
/// When no previous value or seed is available the fallback is empty, so the run
/// proceeds without a coalescent prior. Make that downgrade explicit instead of the
/// misleading "keeping previous Tc" the previous message logged even when there was
/// no previous value to keep.
fn report_coalescent_failure(what: &str, error: &Report, has_fallback: bool) {
  if has_fallback {
    warn!("{what} failed: {error}; retaining the previous coalescent prior for this round");
  } else {
    warn!("{what} failed: {error}; proceeding with no coalescent prior this run");
  }
}

struct PartitionInitResult {
  partitions: PartitionTimetreeAllVec,
  gtr: GTR,
  model_name: GtrModelName,
}

fn initialize_partitions_from_params(
  params: &TimetreeParams,
  graph: &GraphTimetree,
  alphabet: Alphabet,
  aln: Option<&[FastaRecord]>,
) -> Result<PartitionInitResult, Report> {
  let model_name = params.model;

  let aln_data = aln.ok_or_else(|| make_report!("Alignment required for marginal reconstruction"))?;
  let created = create_marginal_partition(graph, 0, alphabet, aln_data, model_name, params.dense)?;

  // Read the GTR from the owning partition (single source of truth), not from a
  // standalone clone. Timetree does not mutate the GTR after creation, so this
  // is behavior-preserving; it keeps the contract consistent with the other
  // pipelines and lets the duplicate `PartitionCreated.gtr` field go away.
  let gtr = match &created.partition {
    MarginalPartition::Sparse(p) => p.gtr().clone(),
    MarginalPartition::Dense(p) => p.gtr().clone(),
  };

  let partition = match created.partition {
    MarginalPartition::Sparse(p) => Arc::new(RwLock::new(PartitionTimetree::Sparse(p))),
    MarginalPartition::Dense(p) => Arc::new(RwLock::new(PartitionTimetree::Dense(p))),
  };

  Ok(PartitionInitResult {
    partitions: vec![partition],
    gtr,
    model_name: created.model_name,
  })
}

fn optimize_branch_lengths_pre_step(
  graph: &GraphTimetree,
  partitions: &[PartitionTimetreeRef],
  no_indels: bool,
) -> Result<(), Report> {
  let old_branch_lengths = save_branch_lengths(graph);

  if no_indels {
    run_optimize_mixed_inner(graph, partitions, BranchOptMethod::BrentSqrt, 0.0, true)
  } else {
    run_optimize_mixed(graph, partitions, BranchOptMethod::BrentSqrt)
  }
  .wrap_err("ML branch-length optimization pre-step failed")?;

  apply_damping(graph, &old_branch_lengths, TIMETREE_PRE_STEP_DAMPING, 0);
  update_marginal(graph, partitions)?;

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::{CoalescentMode, coalescent_mode};
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::disabled(       None,       false, false, CoalescentMode::Disabled)]
  #[case::fixed(          Some(0.25), false, false, CoalescentMode::Fixed(0.25))]
  #[case::opt_default(    None,       true,  false, CoalescentMode::Constant { seed: None })]
  #[case::opt_configured( Some(0.25), true,  false, CoalescentMode::Constant { seed: Some(0.25) })]
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
}
