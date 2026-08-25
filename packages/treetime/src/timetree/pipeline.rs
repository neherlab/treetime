use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::clock::clock_filter::clock_filter_inplace;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
use crate::clock::date_constraints::load_date_constraints;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
use crate::clock::reroot::RerootParams;
use crate::coalescent::coalescent::CoalescentModel;
use crate::coalescent::lineage_counts::compute_lineage_counts;
use crate::coalescent::population_size::effective_population_size;
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
use crate::timetree::inference::runner::{commit_clock_branch_lengths, run_timetree};
use crate::timetree::optimization::clock_filter::{apply_outlier_bad_branches, report_bad_branches};
use crate::timetree::optimization::reroot::reroot_tree;
use crate::timetree::params::{TimeMarginalMode, build_covariation_clock_params, compute_effective_time_marginal};
use crate::timetree::refinement::{Refinement, RefinementOptions, TopologyRefinement};
use crate::timetree::utils::{initialize_clock_totals_from_time_distributions, initialize_node_divergences};
use eyre::{Report, WrapErr};
use log::{debug, info};
use ndarray::array;
use parking_lot::RwLock;
use serde::Serialize;
use std::io::Write;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_io::dates_csv::DatesMap;
use treetime_io::fasta::FastaRecord;
use treetime_utils::make_report;
use treetime_utils::sync::random::get_random_number_generator;

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
  /// Generations per year, used to report the coalescent effective population size
  /// `N_e = Tc * gen_per_year`. Reporting only; does not enter the inference.
  pub gen_per_year: f64,
  pub n_branches_posterior: Option<usize>,
  pub time_marginal: TimeMarginalMode,
  pub confidence: bool,
  pub reconstruct_tip_states: bool,
  pub report_ambiguous: bool,
  pub zero_based: bool,
  /// Seed for the polytomy-resolution sampler. `None` draws one from entropy and logs it,
  /// so a run can be reproduced after the fact.
  pub seed: Option<u64>,
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
  }

  // Initial time tree
  run_timetree(&mut input.graph, &partitions, &clock_model, None, params.no_indels)?;

  // set up coalescent parameters and inference mode
  let skyline_params = SkylineParams {
    n_points: params.n_skyline,
    stiffness: params.skyline_stiffness,
    ..SkylineParams::default()
  };

  if params.n_branches_posterior.is_some() {
    return make_error!("--n-branches-posterior is not yet implemented");
  }
  let coalescent = coalescent_mode(params.coalescent, params.coalescent_opt, params.coalescent_skyline);

  // k(t) frozen for the whole run, read from the first timetree: the earliest point where node
  // times give lineage counts, and before the optimization loop starts editing the topology.
  // Held fixed because it is the prior's input; the statistic role keeps reading the live tree.
  // See `CoalescentModel`.
  let lineage_counts = compute_lineage_counts(&input.graph).wrap_err("Failed to compute coalescent lineage counts")?;

  // Estimate the coalescent Tc (constant or skyline) from the rerooted tree, then re-infer node
  // times under that prior. Estimated after the reroot because the lineage counts it reads are a
  // property of the rooted tree the optimization loop starts from. There is no need to start
  // from a constant Tc and switch to the skyline later: both are cheap analytic solves.
  let mut coalescent_tc = coalescent_timescale(coalescent, &input.graph, &skyline_params)?;

  // Whether that timescale is also a prior on node times. It is estimated either way, because
  // polytomy resolution samples mergers at the per-branch coalescent rate and needs one even for
  // a run that asked for no coalescent.
  let prior_wanted = coalescent != CoalescentMode::Disabled;

  if prior_wanted {
    let prior = CoalescentModel::new(&lineage_counts, &coalescent_tc.distribution)?;
    run_timetree(
      &mut input.graph,
      &partitions,
      &clock_model,
      Some(&prior),
      params.no_indels,
    )?;
  }
  // at this stage we have a consistent coalescent model and timed tree. Subsequence steps are refinement and post-processing.

  // Seed the clock-constrained lengths the loop's first marginal reconstruction propagates along.
  // Undamped: nothing has been committed yet, so there is nothing to blend with.
  commit_clock_branch_lengths(&input.graph, clock_model.clock_rate(), 1.0);

  progress.check_cancelled()?;
  progress.report("Optimization", 0.3, "");
  info!("### TreeTime: Optimisation rounds");
  let mut optimizer = TimetreeOptimizer::new(params.max_iter, false);
  if let Some(writer) = tracelog {
    optimizer = optimizer.with_tracelog(writer)?;
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

  // One generator for the whole loop: polytomy resolution samples a coalescent history per
  // round, and re-seeding per round would correlate them.
  let seed = params.seed.unwrap_or_else(rand::random);
  if params.resolve_polytomies {
    info!("Polytomy resolution is stochastic; seed {seed} (pass --seed to reproduce this run)");
  }
  let mut rng = get_random_number_generator(Some(seed));

  // OPTIMIZATION LOOP
  while let Some(IterationContext { i }) = optimizer.next_iter() {
    progress.check_cancelled()?;
    let iter_fraction = 0.3 + 0.5 * (i as f64 / max_iter as f64);
    progress.report(
      "Optimization",
      iter_fraction,
      &format!("iteration {}/{max_iter}", i + 1),
    );

    // Tc is re-estimated from the live tree -- that is the statistic role, and it is what carries
    // the run's ability to adapt. The lineage counts behind the model stay those of the tree the
    // loop started from.
    if coalescent.is_optimized() {
      coalescent_tc = coalescent_timescale(coalescent, &input.graph, &skyline_params)?;
    }
    let coalescent_model = CoalescentModel::new(&lineage_counts, &coalescent_tc.distribution)?;
    // Preserve every k(t) and Tc(t) discontinuity for event-sampler boundaries.
    let merger_rate = coalescent_model.branch_merger_rate_schedule(&coalescent_tc.schedule)?;

    let outcome = Refinement {
      graph: &mut input.graph,
      partitions: &partitions,
      clock_model: &mut clock_model,
      clock_params: reroot_clock_params,
      branch_params: &branch_params,
      merger_rate: &merger_rate,
      prior: prior_wanted.then_some(&coalescent_model),
      rng: &mut rng,
      options: &refinement_options,
    }
    .run()
    .wrap_err_with(|| format!("When running round {i}"))?;

    optimizer
      .record(
        outcome.sequence_changes,
        outcome.topology.resolved_nodes(),
        outcome.time_change,
        &input.graph,
        &partitions,
        prior_wanted.then_some(&coalescent_tc.distribution),
      )
      .wrap_err("Failed to record convergence metrics")
      .wrap_err_with(|| format!("When running round {i}"))?;
  }

  // Report the effective population size implied by the converged Tc, matching v0's screen output
  // for the constant/opt/skyline coalescent. Fixed and disabled modes report no coalescent, as in
  // v0. N_e = Tc * gen_per_year is a reporting quantity and does not affect the inference.
  if coalescent.is_optimized() {
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

  progress.check_cancelled()?;
  progress.report("Postprocessing", 0.85, "");
  info!("### TreeTime: postprocessing");

  let rate_std = if params.confidence {
    determine_rate_std(params.clock_std_dev, params.covariation, &clock_model)?
  } else {
    None
  };

  // The prior the post-processing passes are inferred under, unchanged across all of them. The
  // susceptibility runs in particular perturb the clock rate, and would otherwise each rebuild
  // k(t) from their own perturbed times, confounding the sensitivity being measured with the
  // prior's reaction to it.
  let final_model = CoalescentModel::new(&lineage_counts, &coalescent_tc.distribution)?;
  let final_prior = prior_wanted.then_some(&final_model);

  if let Some(rate_std) = rate_std {
    info!("### Rate susceptibility analysis (rate_std={rate_std:.6e})");
    compute_rate_susceptibility(
      &mut input.graph,
      &partitions,
      &clock_model,
      final_prior,
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
      final_prior,
      params.no_indels,
    )
    .wrap_err("Final timetree inference failed")?;

    // Undamped: this reconstruction reports the final tree, so it runs on the lengths these
    // final times imply rather than on a blend with the loop's last round.
    commit_clock_branch_lengths(&input.graph, clock_model.clock_rate(), 1.0);

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
  /// Optimize a constant Tc each round via the exact analytic solve.
  Constant,
  /// Optimize a piecewise-constant skyline Tc(t) each round.
  Skyline,
}

impl CoalescentMode {
  /// Whether Tc is re-estimated from the tree, as opposed to fixed or disabled.
  fn is_optimized(self) -> bool {
    matches!(self, CoalescentMode::Constant | CoalescentMode::Skyline)
  }
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

/// Estimates the coalescent Tc prior for the current tree under `mode`.
///
/// Both the constant and skyline solves are exact and cheap, so this runs every
/// round. They fail only on a tree that is degenerate for the coalescent (no time
/// span or no mergers); such a failure is propagated so the run stops with a clear
/// message rather than silently substituting an invented timescale.
fn estimate_coalescent_tc(
  mode: CoalescentMode,
  graph: &GraphTimetree,
  skyline_params: &SkylineParams,
) -> Result<Option<CoalescentTimescale>, Report> {
  // Both optimizing modes are the same solve: a constant Tc is just a one-segment
  // skyline. They differ only in the number of segments.
  let n_points = match mode {
    CoalescentMode::Disabled => return Ok(None),
    CoalescentMode::Fixed(tc) => return Ok(Some(CoalescentTimescale::constant(tc))),
    CoalescentMode::Constant => 1,
    CoalescentMode::Skyline => skyline_params.n_points,
  };
  let result = optimize_skyline(
    graph,
    &SkylineParams {
      n_points,
      ..skyline_params.clone()
    },
  )?;
  Ok(Some(CoalescentTimescale {
    distribution: result.tc_distribution,
    schedule: result.tc_schedule,
  }))
}

/// The coalescent timescale a run works at, whether or not it asked for a coalescent prior.
///
/// Polytomy resolution samples mergers at the per-branch coalescent rate, so a Tc is needed even
/// under [`CoalescentMode::Disabled`]. There it comes from the same one-segment analytic solve
/// `--coalescent-opt` uses, rather than from v0's dummy rate, which is calibrated per polytomy to
/// the very time window the sampled history has to fit into.
fn coalescent_timescale(
  mode: CoalescentMode,
  graph: &GraphTimetree,
  skyline_params: &SkylineParams,
) -> Result<CoalescentTimescale, Report> {
  let mode = match mode {
    CoalescentMode::Disabled => CoalescentMode::Constant,
    mode => mode,
  };
  estimate_coalescent_tc(mode, graph, skyline_params)
    .wrap_err("Failed to estimate the coalescent timescale")?
    .ok_or_else(|| make_report!("A coalescent Tc is required, but {mode:?} yielded none"))
}

struct CoalescentTimescale {
  distribution: Distribution,
  schedule: PiecewiseConstantFn,
}

impl CoalescentTimescale {
  fn constant(tc: f64) -> Self {
    Self {
      distribution: Distribution::constant(tc),
      schedule: PiecewiseConstantFn::new(array![], array![tc]),
    }
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
}
