use crate::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::clock::clock_filter::clock_filter_inplace;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_output::write_clock_model;
use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
use crate::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::clock::reroot::RerootParams;
use crate::commands::shared::args::BranchLengthMode;
use crate::commands::timetree::args::{TimeMarginalMode, TreetimeTimetreeArgs};
use crate::commands::timetree::coalescent::optimize_tc::optimize_tc;
use crate::commands::timetree::coalescent::skyline::{SkylineParams, optimize_skyline};
use crate::commands::timetree::convergence::metrics::{IterationContext, TimetreeOptimizer};
use crate::commands::timetree::inference::runner::run_timetree;
use crate::commands::timetree::initialization::{InputData, initialize_partitions, load_input_data};
use crate::commands::timetree::optimization::clock_filter::{apply_outlier_bad_branches, report_bad_branches};
use crate::commands::timetree::optimization::reroot::reroot_tree;
use crate::commands::timetree::output::auspice::write_auspice_json;
use crate::commands::timetree::output::confidence::{
  NodeConfidenceInterval, compute_rate_susceptibility, determine_rate_std, extract_confidence_intervals,
  write_confidence_intervals,
};
use crate::commands::timetree::refinement::run_refinement_iteration;
use crate::commands::timetree::utils::initialize_clock_totals_from_time_distributions;
use crate::optimize::args::BranchOptMethod;
use crate::optimize::iteration::{apply_damping, save_branch_lengths};
use crate::optimize::optimize_unified::{run_optimize_mixed, run_optimize_mixed_inner};
use crate::representation::partition::timetree::GraphTimetree;
use crate::representation::partition::traits::PartitionTimetreeAll;
use crate::representation::payload::ancestral::annotate_branch_mutations;
use crate::representation::payload::timetree::EdgeTimetree;
use crate::representation::payload::timetree::NodeTimetree;
use crate::seq::alignment::get_common_length;
use crate::{make_error, make_report};
use eyre::{Report, WrapErr};
use log::{debug, info, warn};
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_io::fasta::FastaRecord;
use treetime_io::graph::write_graph_files;

/// v0 default damping for the timetree pre-step branch-length optimization (treeanc.py:1298).
const TIMETREE_PRE_STEP_DAMPING: f64 = 0.75;

pub fn run_timetree_estimation(args: &TreetimeTimetreeArgs) -> Result<(), Report> {
  info!("# TreeTime Timetree Estimation");
  debug!(
    "Branch length mode: {:?}, Keep root: {}",
    args.branch_length_mode, args.keep_root
  );

  // Compute effective time_marginal, accounting for --confidence flag.
  //
  // v0 (wrappers.py:453-469):
  //   - --confidence with --covariation or --clock-std-dev: set time_marginal='confidence-only'
  //     (equivalent to v1's OnlyFinal -- runs a final marginal pass for CI estimation)
  //   - --confidence without prerequisites: warn and disable confidence
  let time_marginal = compute_effective_time_marginal(args);

  info!("## Loading input data");
  let InputData {
    mut graph,
    alphabet,
    aln,
  } = load_input_data(args)?;
  info!(
    "Input loaded: {} nodes, {} edges",
    graph.get_nodes().len(),
    graph.get_edges().len()
  );

  // Covariation-aware clock regression: weighted least squares with phylogenetic
  // variance structure (Neher 2018, "Efficient estimation of evolutionary rates
  // by covariance aware regression").
  //
  // v0 (clock_tree.py:277-285):
  //   branch_variance = (max(0, clock_length) + tip_slack^2 * om) * om  [leaves]
  //   branch_variance = max(0, clock_length) * om                       [internal]
  // where om = one_mutation = 1/seq_len, tip_slack = OVER_DISPERSION = 10
  //
  // Without covariation (default): OLS with unit variance for leaves, zero for internal.
  // This ignores phylogenetic correlation between tips, producing unreliable CIs.
  //
  // v0 uses default (non-covariation) params for the initial regression and clock
  // filter, switching to covariation params only during iteration-loop reroots
  // (treetime.py:488,491 vs 643-644). We follow the same pattern.
  let covariation_clock_params = build_covariation_clock_params(args, aln.as_deref())?;
  let branch_params = BranchPointOptimizationParams::default();

  // Initial regression: always non-covariation (matching v0 treetime.py:488-491)
  let reroot_params = RerootParams {
    force_positive_rate: !args.allow_negative_rate,
    ..RerootParams::default()
  };
  let mut clock_model = estimate_clock_model_with_reroot_policy(
    &mut graph,
    &ClockParams::default(),
    args.clock_rate,
    args.keep_root,
    &branch_params,
    &reroot_params,
    None,
  )
  .wrap_err("Failed to infer clock model")?
  .clock_model;

  let partitions = match args.branch_length_mode {
    BranchLengthMode::Input => {
      info!("Branch length mode: Input - using tree branch lengths");
      vec![]
    },
    BranchLengthMode::Marginal => {
      info!("Branch length mode: Marginal - initializing partitions from alignment");
      initialize_partitions(args, &graph, alphabet, aln.as_deref())?
    },
  };

  // ML branch-length optimization pre-step (v0: optimize_tree(max_iter=1))
  //
  // v0 calls optimize_tree(infer_gtr=True, max_iter=1) before rerooting,
  // then optimize_tree(max_iter=1) again after rerooting + clock filter.
  // Both calls run marginal reconstruction followed by one round of per-edge
  // Brent optimization, seeding time inference with ML-optimized branch
  // lengths rather than raw input values.
  if let Some(aln) = aln.as_deref() {
    if args.branch_length_mode == BranchLengthMode::Marginal && !partitions.is_empty() {
      info!("### ML branch-length optimization (pre-reroot)");
      initialize_marginal(&graph, &partitions, aln)?;
      optimize_branch_lengths_pre_step(&graph, &partitions, args.no_indels)
        .wrap_err("ML branch-length optimization (pre-reroot) failed")?;
    }
  }

  // First reroot: use non-covariation params (pre-filter, matching v0)
  if !args.keep_root {
    info!("First reroot (pre-ancestral)");
    clock_model = reroot_tree(
      &mut graph,
      &partitions,
      &ClockParams::default(),
      args.clock_rate,
      &branch_params,
      !args.allow_negative_rate,
    )
    .wrap_err("Failed to reroot tree (pre-ancestral)")?;
  }

  if args.clock_filter > 0.0 {
    let result = clock_filter_inplace(&graph, &clock_model, args.clock_filter);
    report_bad_branches(&graph, &clock_model, result.iqd);
    apply_outlier_bad_branches(&graph);
  }

  if let Some(aln) = aln.as_deref() {
    match args.branch_length_mode {
      BranchLengthMode::Input => {
        info!("Using input branch lengths for timetree inference");
      },
      BranchLengthMode::Marginal => {
        info!("### ML branch-length optimization (post-reroot)");
        update_marginal(&graph, &partitions)?;
        optimize_branch_lengths_pre_step(&graph, &partitions, args.no_indels)
          .wrap_err("ML branch-length optimization (post-reroot) failed")?;
      },
    }
  }

  info!("### TreeTime: initial round");
  info!("### Initializing node times from date constraints");
  initialize_clock_totals_from_time_distributions(&graph)?;

  // Skyline coalescent parameters (constructed once, reused for initial and final optimization)
  let skyline_params = SkylineParams {
    n_points: args.n_skyline,
    ..SkylineParams::default()
  };

  if args.n_branches_posterior.is_some() {
    return make_error!("--n-branches-posterior is not yet implemented");
  }

  // First pass without coalescent: establish internal node time distributions via
  // backward+forward pass. Coalescent contributions read node times that only exist
  // after the backward pass writes them.
  run_timetree(&mut graph, &partitions, &clock_model, None, args.no_indels)?;

  // Build coalescent Tc using established node times.
  // For skyline mode: use constant Tc during iterations, switch to skyline after convergence.
  let mut coalescent_tc: Option<Distribution> = if args.coalescent_skyline {
    // Starting point for Brent's method bounded search; insensitive to initial value
    let initial_tc = 1.0;
    info!("### Optimizing constant coalescent Tc (pre-loop, skyline deferred)");
    match optimize_tc(&graph, initial_tc) {
      Ok(result) if result.success => {
        info!(
          "Pre-loop Tc = {:.6e} (likelihood = {:.4})",
          result.tc, result.likelihood
        );
        Some(Distribution::constant(result.tc))
      },
      Ok(_) => {
        warn!("Pre-loop Tc optimization did not converge, using Tc = {initial_tc:.6e}");
        Some(Distribution::constant(initial_tc))
      },
      Err(e) => {
        warn!("Pre-loop Tc optimization failed: {e}, using Tc = {initial_tc:.6e}");
        Some(Distribution::constant(initial_tc))
      },
    }
  } else {
    args.coalescent.map(Distribution::constant)
  };

  // Second pass with coalescent prior
  if coalescent_tc.is_some() {
    run_timetree(
      &mut graph,
      &partitions,
      &clock_model,
      coalescent_tc.as_ref(),
      args.no_indels,
    )?;
  }

  // Post-filter reroots: use covariation params when --covariation is enabled,
  // matching v0's pattern where covariation kicks in from the iteration loop onward
  // (treetime.py:643-644).
  let default_clock_params = ClockParams::default();
  let reroot_clock_params = covariation_clock_params.as_ref().unwrap_or(&default_clock_params);

  if !args.keep_root {
    info!("Reroot (post-ancestral)");
    clock_model = reroot_tree(
      &mut graph,
      &partitions,
      reroot_clock_params,
      args.clock_rate,
      &branch_params,
      !args.allow_negative_rate,
    )
    .wrap_err("Failed to reroot tree (post-ancestral)")?;
  }

  info!("### TreeTime: Optimisation rounds");
  let mut optimizer = TimetreeOptimizer::new(args.max_iter, args.coalescent_skyline, args.tracelog.clone())?;
  while let Some(IterationContext { i }) = optimizer.next_iter() {
    // Re-optimize constant Tc each iteration. Skip i < 2 because v1 runs a pre-loop
    // Tc optimization, making the first in-loop re-optimization redundant.
    // In skyline mode, constant Tc is used during loop iterations; the full skyline
    // fit is deferred to post-convergence.
    if (args.coalescent_opt || args.coalescent_skyline) && i >= 2 {
      // For constant distributions, max_value() returns the Tc amplitude
      let initial_tc = coalescent_tc.as_ref().map_or(1.0, |d| d.max_value());
      match optimize_tc(&graph, initial_tc) {
        Ok(result) => {
          if result.success {
            coalescent_tc = Some(Distribution::constant(result.tc));
            info!(
              "Optimized Tc = {:.6e} (likelihood = {:.4})",
              result.tc, result.likelihood
            );
          } else {
            warn!("Tc optimization did not converge, keeping Tc = {initial_tc:.6e}");
          }
        },
        Err(e) => {
          warn!("Tc optimization failed: {e}, keeping Tc = {initial_tc:.6e}");
        },
      }
    }

    let (n_diff, n_resolved) = run_refinement_iteration(
      args,
      &mut graph,
      &partitions,
      &mut clock_model,
      reroot_clock_params,
      &branch_params,
      coalescent_tc.as_ref(),
    )
    .wrap_err_with(|| format!("When running round {i}"))?;

    optimizer
      .record(n_diff, n_resolved, &graph, &partitions, coalescent_tc.as_ref())
      .wrap_err("Failed to record convergence metrics")
      .wrap_err_with(|| format!("When running round {i}"))?;
  }

  // Re-optimize skyline with stabilized node times (matching v0 behavior where skyline
  // optimization happens only on the final iteration after times have converged)
  if args.coalescent_skyline {
    info!("### Re-optimizing skyline coalescent with stabilized node times");
    let skyline_result =
      optimize_skyline(&graph, &skyline_params).wrap_err("Failed to re-optimize skyline coalescent model")?;
    info!(
      "Skyline re-optimization completed: log_likelihood={:.4}",
      skyline_result.log_likelihood
    );
    coalescent_tc = Some(skyline_result.tc_distribution);

    // Run final timetree pass with optimized skyline. OnlyFinal defers this to the
    // final marginal pass below; Always and Never take this path.
    if time_marginal != TimeMarginalMode::OnlyFinal {
      run_timetree(
        &mut graph,
        &partitions,
        &clock_model,
        coalescent_tc.as_ref(),
        args.no_indels,
      )
      .wrap_err("Final timetree pass with optimized skyline failed")?;

      if !partitions.is_empty() {
        update_marginal(&graph, &partitions)?;
      }
    }
  }

  info!("### TreeTime: postprocessing");

  // Rate susceptibility analysis: re-run timetree at rate +/- rate_std.
  // Runs BEFORE the final marginal pass (matching v0 treetime.py:380-388).
  //
  // v0 gating (wrappers.py:453-469): rate susceptibility requires --confidence
  // AND either --clock-std-dev or --covariation. Without --confidence, rate
  // susceptibility never runs even if rate_std is available.
  //
  // v1's run_timetree always produces marginal distributions. Rate susceptibility
  // runs three passes (upper, lower, central rate); the central-rate pass restores
  // the graph to its pre-call state.
  let rate_std = if args.confidence {
    determine_rate_std(args.clock_std_dev, args.covariation, &clock_model)?
  } else {
    None
  };

  if let Some(rate_std) = rate_std {
    let current_rate = clock_model.clock_rate();
    if current_rate <= 0.0 {
      warn!("Clock rate is non-positive ({current_rate:.6e}), skipping rate susceptibility");
    } else {
      info!("### Rate susceptibility analysis (rate_std={rate_std:.6e})");
      compute_rate_susceptibility(
        &mut graph,
        &partitions,
        &clock_model,
        coalescent_tc.as_ref(),
        rate_std,
        args.no_indels,
      )
      .wrap_err("Rate susceptibility analysis failed")?;
    }
  }

  // Final marginal pass for confidence interval estimation
  if time_marginal == TimeMarginalMode::OnlyFinal {
    info!("### Final round: marginal reconstruction for confidence intervals");
    run_timetree(
      &mut graph,
      &partitions,
      &clock_model,
      coalescent_tc.as_ref(),
      args.no_indels,
    )
    .wrap_err("Final timetree inference failed")?;

    if !partitions.is_empty() {
      update_marginal(&graph, &partitions)?;
    }
  }

  // Extract CI when any CI data is available.
  //
  // v0 behavior: marginal distributions provide mutation-stochasticity CIs via HPD
  // regions. OnlyFinal and Always both produce these distributions (v1 always uses
  // marginal inference, so both modes have HPD data). Rate susceptibility adds
  // clock-rate-uncertainty CIs independently.
  let confidence_intervals =
    if matches!(time_marginal, TimeMarginalMode::OnlyFinal | TimeMarginalMode::Always) || rate_std.is_some() {
      let intervals = extract_confidence_intervals(&graph);
      let ci_path = args.outdir.join("confidence_intervals.tsv");
      write_confidence_intervals(&intervals, &ci_path).wrap_err("Failed to write confidence intervals")?;
      info!("Wrote confidence intervals to {ci_path}", ci_path = ci_path.display());
      Some(intervals)
    } else {
      None
    };

  info!("### TreeTime: writing outputs");
  write_outputs(args, &graph, &partitions, &clock_model, confidence_intervals.as_deref())?;

  Ok(())
}

/// Compute effective time_marginal mode, promoting Never to OnlyFinal when
/// --confidence is set with rate uncertainty prerequisites.
///
/// v0 (wrappers.py:478):
///   time_marginal = 'confidence-only' if (calc_confidence and time_marginal == 'never') else time_marginal
fn compute_effective_time_marginal(args: &TreetimeTimetreeArgs) -> TimeMarginalMode {
  if args.confidence && args.time_marginal == TimeMarginalMode::Never {
    if args.clock_std_dev.is_some() || args.covariation {
      info!("--confidence: promoting time-marginal from never to only-final for CI estimation");
      TimeMarginalMode::OnlyFinal
    } else {
      warn!(
        "Cannot estimate confidence intervals without clock rate uncertainty. \
         Specify --clock-std-dev or rerun with --covariation. \
         Proceeding without confidence estimation."
      );
      TimeMarginalMode::Never
    }
  } else {
    args.time_marginal
  }
}

/// Build covariation-aware ClockParams when --covariation is enabled.
///
/// v0 (clock_tree.py:277-285):
///   branch_variance = (max(0, clock_length) + tip_slack^2 * om) * om   [leaves]
///   branch_variance = max(0, clock_length) * om                         [internal]
/// where om = 1/seq_len, tip_slack = OVER_DISPERSION = 10 (config.py:8)
///
/// In v1's ClockParams model:
///   branch_variance = variance_factor * edge_len + variance_offset        [all]
///   branch_variance += variance_offset_leaf                               [leaves only]
///
/// Mapping:
///   variance_factor = 1.0 / seq_len   (= om, the one_mutation value)
///   variance_offset = 0.0
///   variance_offset_leaf = tip_slack^2 / seq_len^2   (= tip_slack^2 * om^2)
fn build_covariation_clock_params(
  args: &TreetimeTimetreeArgs,
  aln: Option<&[FastaRecord]>,
) -> Result<Option<ClockParams>, Report> {
  if !args.covariation {
    return Ok(None);
  }

  let seq_len = if let Some(aln_data) = aln {
    get_common_length(aln_data)? as f64
  } else {
    args
      .sequence_length
      .ok_or_else(|| make_report!("--sequence-length required for --covariation without alignment"))? as f64
  };

  // v0 default: OVER_DISPERSION = 10 (config.py:8), overridable via --tip-slack
  let tip_slack = args.tip_slack.unwrap_or(10.0);

  info!("Covariation-aware clock regression: seq_len={seq_len}, tip_slack={tip_slack}");

  Ok(Some(ClockParams {
    variance_factor: 1.0 / seq_len,
    variance_offset: 0.0,
    variance_offset_leaf: tip_slack * tip_slack / (seq_len * seq_len),
  }))
}

/// Run one pass of ML branch-length optimization on timetree partitions.
///
/// v0: `optimize_tree(max_iter=1)` calls `optimize_tree_marginal(max_iter=1,
/// damping=0.75)`, which does per-edge Brent optimization in sqrt(t) space
/// with damped updates followed by marginal reconstruction. At iteration 0
/// with damping=0.75, the update blends 25% optimized + 75% old branch
/// lengths (`damping_factor = 0.75^1 = 0.75`).
///
/// Requires marginal profiles to be populated (via `initialize_marginal` or
/// `update_marginal`) before calling.
fn optimize_branch_lengths_pre_step(
  graph: &GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  no_indels: bool,
) -> Result<(), Report> {
  let old_branch_lengths = save_branch_lengths(graph);

  if no_indels {
    run_optimize_mixed_inner(graph, partitions, BranchOptMethod::BrentSqrt, 0.0, true)
  } else {
    run_optimize_mixed(graph, partitions, BranchOptMethod::BrentSqrt)
  }
  .wrap_err("ML branch-length optimization pre-step failed")?;

  // Blend optimized branch lengths with old values (iteration=0):
  // bl = bl_optimized * (1 - 0.75^1) + bl_old * 0.75^1
  //    = bl_optimized * 0.25 + bl_old * 0.75
  apply_damping(graph, &old_branch_lengths, TIMETREE_PRE_STEP_DAMPING, 0);

  // Re-run marginal reconstruction with damped branch lengths
  update_marginal(graph, partitions)?;

  Ok(())
}

fn write_outputs(
  args: &TreetimeTimetreeArgs,
  graph: &GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  clock_model: &ClockModel,
  confidence_intervals: Option<&[NodeConfidenceInterval]>,
) -> Result<(), Report> {
  // Populate the `mutations` comment on every branch from the current partition
  // state before tree serialization. v0 emits `[&mutations="..."]` on every
  // node in the Nexus tree (CLI_io.py:167-199); v1 now matches that via the
  // `NodeAncestral.mutations` field which `NodeTimetree.nwk_comments()`
  // inherits from its base payload. Partitions must have completed marginal
  // reconstruction before this call so that `edge_subs()` reads the final
  // state rather than stale Fitch parsimony data.
  if !partitions.is_empty() {
    annotate_branch_mutations(graph, partitions).wrap_err("Failed to annotate branch mutations for tree output")?;
  }

  write_graph_files(&args.outdir, "timetree", graph).wrap_err("Failed to write tree output")?;

  write_clock_model(clock_model, &args.outdir.join("timetree"))?;

  write_auspice_json(graph, confidence_intervals, &args.outdir)?;

  if args.plot_rtt.is_some() {
    return make_error!("--plot-rtt is not yet implemented");
  }

  if args.plot_tree.is_some() {
    return make_error!("--plot-tree is not yet implemented");
  }

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  fn args_with(
    time_marginal: TimeMarginalMode,
    confidence: bool,
    covariation: bool,
    clock_std_dev: Option<f64>,
  ) -> TreetimeTimetreeArgs {
    TreetimeTimetreeArgs {
      clock_std_dev,
      time_marginal,
      confidence,
      covariation,
      ..TreetimeTimetreeArgs::default()
    }
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::never_no_confidence(            (TimeMarginalMode::Never,     false, false, None),      TimeMarginalMode::Never)]
  #[case::never_confidence_covariation(   (TimeMarginalMode::Never,     true,  true,  None),      TimeMarginalMode::OnlyFinal)]
  #[case::never_confidence_clock_std(     (TimeMarginalMode::Never,     true,  false, Some(0.1)), TimeMarginalMode::OnlyFinal)]
  #[case::never_confidence_no_prereqs(    (TimeMarginalMode::Never,     true,  false, None),      TimeMarginalMode::Never)]
  #[case::always_no_confidence(           (TimeMarginalMode::Always,    false, false, None),      TimeMarginalMode::Always)]
  #[case::always_confidence_covariation(  (TimeMarginalMode::Always,    true,  true,  None),      TimeMarginalMode::Always)]
  #[case::always_confidence_clock_std(    (TimeMarginalMode::Always,    true,  false, Some(0.1)), TimeMarginalMode::Always)]
  #[case::only_final_no_confidence(       (TimeMarginalMode::OnlyFinal, false, false, None),      TimeMarginalMode::OnlyFinal)]
  #[case::only_final_confidence(          (TimeMarginalMode::OnlyFinal, true,  true,  None),      TimeMarginalMode::OnlyFinal)]
  #[trace]
  fn test_compute_effective_time_marginal(
    #[case] (mode, confidence, covariation, clock_std_dev): (TimeMarginalMode, bool, bool, Option<f64>),
    #[case] expected: TimeMarginalMode,
  ) {
    let args = args_with(mode, confidence, covariation, clock_std_dev);
    assert_eq!(expected, compute_effective_time_marginal(&args));
  }
}
