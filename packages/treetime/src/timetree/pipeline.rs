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
use crate::commands::timetree::output::coalescent::{
  CoalescentBand, CoalescentInputs, CoalescentOutput, CoalescentOutputMode, CoalescentSolve,
};
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::make_error;
use crate::optimize::dispatch::{run_optimize_mixed, run_optimize_mixed_inner};
use crate::optimize::iteration::{apply_damping, save_branch_lengths};
use crate::optimize::params::{BranchLengthMode, BranchOptMethod};
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::timetree::partition::{
  GraphTimetree, PartitionTimetree, PartitionTimetreeAllVec, PartitionTimetreeRef,
};
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
use ndarray::{Array1, array};
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
  pub skyline_n_points: usize,
  /// Skyline smoothing stiffness (penalizes changes in 1/Tc; units of time^2).
  pub skyline_stiffness: f64,
  /// Confidence level (standard deviations) for inferred coalescent Tc bands.
  pub coalescent_confidence: f64,
  /// Generations per year, used to report the coalescent effective population size
  /// `N_e = Tc * gen_per_year`. Reporting only; does not enter the inference.
  pub gen_per_year: f64,
  pub n_branches_posterior: Option<usize>,
  pub time_marginal: TimeMarginalMode,
  pub confidence: bool,
  pub include_leaves: bool,
  pub impute_missing_data: bool,
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
  /// The inferred coalescent time scale, or `None` when the run asked for no coalescent.
  #[serde(skip)]
  pub coalescent: Option<CoalescentOutput>,
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
    n_points: params.skyline_n_points,
    stiffness: params.skyline_stiffness,
    n_std: params.coalescent_confidence,
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

  // Report the effective population size implied by the converged Tc on screen, for every mode that
  // writes a coalescent output file so the two channels agree. This extends v0's constant/opt/skyline
  // screen output to a fixed Tc; only the disabled mode reports nothing. N_e = Tc * gen_per_year is a
  // reporting quantity and does not affect the inference.
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

  let coalescent_output = build_coalescent_output(coalescent, &coalescent_tc, params.gen_per_year, &skyline_params)?;

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

  /// Serialization tag for the coalescent output document, or `None` only when the run writes no
  /// coalescent output (`Disabled`). Every set coalescent is written -- a fixed user Tc as well as
  /// an inferred constant or skyline. Emitting a fixed Tc diverges from v0, which writes nothing for
  /// it, by explicit decision (see `kb/decisions/coalescent-output-schema.md`). The stderr `N_e`
  /// report uses the same gate, so the file and screen agree. Sole mapping from the inference mode
  /// to the output tag, so a new mode is one edit here.
  fn output_mode(self) -> Option<CoalescentOutputMode> {
    match self {
      CoalescentMode::Disabled => None,
      CoalescentMode::Fixed(_) => Some(CoalescentOutputMode::Fixed),
      CoalescentMode::Constant => Some(CoalescentOutputMode::Constant),
      CoalescentMode::Skyline => Some(CoalescentOutputMode::Skyline),
    }
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
    CoalescentMode::Fixed(tc) => return fixed_timescale(tc, graph).map(Some),
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

/// Builds the timescale for a fixed, user-supplied Tc: a constant schedule plus a one-segment report
/// spanning the tree's full time range.
///
/// The span is the range of the coalescent breakpoints -- the same node-time range the optimized
/// modes report -- so a fixed Tc reports the whole tree as its single segment. No solve runs (a fixed
/// Tc is the escape hatch for trees the optimizer rejects as degenerate), so the report carries no
/// confidence band and no likelihood.
fn fixed_timescale(tc: f64, graph: &GraphTimetree) -> Result<CoalescentTimescale, Report> {
  let lineage_counts = compute_lineage_counts(graph).wrap_err("Failed to compute coalescent lineage counts")?;
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
  /// Per-segment reporting data (boundaries and confidence band) from the analytic solve. `None`
  /// for a fixed user Tc, which has no solve and no band.
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

/// Per-segment reporting data carried out of the skyline/constant analytic solve, or synthesized for
/// a fixed Tc, for the output document. Holds the segment boundaries and, for an inferred Tc, the
/// confidence band and likelihood; the Tc values themselves come from
/// [`CoalescentTimescale::schedule`].
struct CoalescentTcReport {
  /// Segment boundaries in numeric date (length `n_segments + 1`, ascending).
  segment_boundaries: Array1<f64>,
  /// Per-segment confidence band, or `None` for a fixed Tc (no solve, no band).
  band: Option<CoalescentReportBand>,
  /// Coalescent log-likelihood at the reported Tc, or `None` for a fixed Tc (not inferred).
  log_likelihood: Option<f64>,
}

/// Per-segment Tc confidence band carried out of the analytic solve. Both bounds are present
/// together, so a partial band is unrepresentable.
struct CoalescentReportBand {
  /// Lower confidence bound per segment.
  lower: Array1<f64>,
  /// Upper confidence bound per segment.
  upper: Array1<f64>,
}

/// Assembles the coalescent output document for the requested mode, or `None` when the run writes no
/// coalescent output (disabled, per [`CoalescentMode::output_mode`]).
///
/// An inferred Tc (constant, skyline) carries the confidence band, likelihood, and per-segment
/// boundaries from the analytic solve. A fixed Tc carries one segment over the tree span with no band
/// and no likelihood. The skyline grid inputs are recorded only for a skyline.
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

  // The skyline grid and stiffness apply only to a skyline; the confidence width only when a band
  // was estimated (a fixed Tc has neither).
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
  use super::{
    CoalescentBand, CoalescentInputs, CoalescentMode, CoalescentOutput, CoalescentOutputMode, CoalescentReportBand,
    CoalescentSolve, CoalescentTcReport, CoalescentTimescale, build_coalescent_output, coalescent_mode,
    estimate_coalescent_tc,
  };
  use crate::clock::date_constraints::load_date_constraints;
  use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
  use crate::partition::timetree::partition::GraphTimetree;
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_distribution::Distribution;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::{o, pretty_assert_array_eq};

  // N_e = T_c * gen_per_year (packages/treetime/src/coalescent/population_size.rs).
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

  // Oracle: the independently tested `CoalescentOutput::new`
  // (packages/treetime/src/commands/timetree/output/__tests__/test_coalescent_output.rs). These
  // tests fix the pipeline-to-output mapping: requested mode, band presence, and (for a fixed Tc)
  // the single-segment span taken from the lineage-count breakpoints.

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
    // A fixed user Tc now writes a single-segment, band-less document spanning the tree, diverging
    // from v0 by explicit decision (kb/decisions/coalescent-output-schema.md). The span is the
    // lineage-count breakpoint range, matching the tree's [2000, 2010] date span.
    let graph = dated_tree()?;
    let params = SkylineParams {
      n_std: N_STD,
      ..SkylineParams::default()
    };
    let timescale = estimate_coalescent_tc(CoalescentMode::Fixed(2.5), &graph, &params)?
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

  fn dated_tree() -> Result<GraphTimetree, Report> {
    // Small dated 3-tip tree spanning [2000, 2010] with two binary mergers, enough for a
    // multi-segment skyline solve.
    let dates: DatesMap = btreemap! {
      o!("root") => Some(DateConstraint::exact(2000.0)),
      o!("x")    => Some(DateConstraint::exact(2005.0)),
      o!("a")    => Some(DateConstraint::exact(2010.0)),
      o!("b")    => Some(DateConstraint::exact(2010.0)),
      o!("c")    => Some(DateConstraint::exact(2010.0)),
    };
    let graph = nwk_read_str("((a:1,b:1)x:1,c:1)root:0;")?;
    load_date_constraints(&dates, &graph)?;
    Ok(graph)
  }

  #[test]
  fn test_pipeline_estimate_coalescent_tc_report_carries_the_skyline_solve() -> Result<(), Report> {
    // The end-to-end pipeline test is disabled on the mass-sizing bug, so the link from the
    // optimizer result to the report struct has no running guard. This exercises exactly that
    // copy: `estimate_coalescent_tc` must carry the solve's own boundaries and band into the
    // `CoalescentTcReport`, unpermuted and unresized. The oracle is a direct `optimize_skyline`
    // call with the same deterministic inputs, so any drop or transpose in the copy shows up.
    let graph = dated_tree()?;
    let params = SkylineParams {
      n_points: 3,
      ..SkylineParams::default()
    };

    let solve = optimize_skyline(&graph, &params)?;
    let timescale = estimate_coalescent_tc(CoalescentMode::Skyline, &graph, &params)?
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
