use crate::commands::ancestral::fitch::get_common_length;
use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::commands::clock::clock_filter::clock_filter_inplace;
use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_output::write_clock_model;
use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::clock::reroot::RerootParams;
use crate::commands::timetree::args::{BranchLengthMode, TimeMarginalMode, TreetimeTimetreeArgs};
use crate::commands::timetree::coalescent::optimize_tc::optimize_tc;
use crate::commands::timetree::coalescent::skyline::{SkylineParams, optimize_skyline};
use crate::commands::timetree::convergence::metrics::{IterationContext, TimetreeOptimizer};
use crate::commands::timetree::inference::runner::run_timetree;
use crate::commands::timetree::initialization::{InputData, initialize_partitions, load_input_data};
use crate::commands::timetree::optimization::clock_filter::{apply_outlier_bad_branches, report_bad_branches};
use crate::commands::timetree::optimization::reroot::reroot_tree;
use crate::commands::timetree::output::confidence::{
  compute_rate_susceptibility, determine_rate_std, extract_confidence_intervals, write_confidence_intervals,
};
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::commands::timetree::refinement::run_refinement_iteration;
use crate::commands::timetree::utils::initialize_clock_totals_from_time_distributions;
use crate::make_error;
use crate::representation::partition::timetree::GraphTimetree;
use crate::representation::payload::timetree::EdgeTimetree;
use crate::representation::payload::timetree::NodeTimetree;
use eyre::{Report, WrapErr};
use log::{debug, info, warn};
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_io::fasta::FastaRecord;
use treetime_io::nex::{NexWriteOptions, nex_write_file};
use treetime_io::nwk::{NwkWriteOptions, nwk_write_file};

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
        info!("Running initial ancestral reconstruction");
        initialize_marginal(&graph, &partitions, aln)?;
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
  run_timetree(&mut graph, &partitions, &clock_model, None)?;

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
    run_timetree(&mut graph, &partitions, &clock_model, coalescent_tc.as_ref())?;
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
      .record(n_diff, n_resolved, &graph, &partitions)
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
      run_timetree(&mut graph, &partitions, &clock_model, coalescent_tc.as_ref())
        .wrap_err("Final timetree pass with optimized skyline failed")?;

      if !partitions.is_empty() {
        update_marginal(&graph, &partitions)?;
      }
    }
  }

  // --- Postprocessing ---

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
      compute_rate_susceptibility(&mut graph, &partitions, &clock_model, coalescent_tc.as_ref(), rate_std)
        .wrap_err("Rate susceptibility analysis failed")?;
    }
  }

  // Final marginal pass for confidence interval estimation
  if time_marginal == TimeMarginalMode::OnlyFinal {
    info!("### Final round: marginal reconstruction for confidence intervals");
    run_timetree(&mut graph, &partitions, &clock_model, coalescent_tc.as_ref())
      .wrap_err("Final timetree inference failed")?;

    if !partitions.is_empty() {
      update_marginal(&graph, &partitions)?;
    }
  }

  // Write CI file when any CI data is available.
  //
  // v0 behavior: marginal distributions provide mutation-stochasticity CIs via HPD
  // regions. OnlyFinal and Always both produce these distributions (v1 always uses
  // marginal inference, so both modes have HPD data). Rate susceptibility adds
  // clock-rate-uncertainty CIs independently.
  if matches!(time_marginal, TimeMarginalMode::OnlyFinal | TimeMarginalMode::Always) || rate_std.is_some() {
    let intervals = extract_confidence_intervals(&graph);
    let ci_path = args.outdir.join("confidence_intervals.tsv");
    write_confidence_intervals(&intervals, &ci_path).wrap_err("Failed to write confidence intervals")?;
    info!("Wrote confidence intervals to {ci_path}", ci_path = ci_path.display());
  }

  info!("### TreeTime: writing outputs");
  write_outputs(args, &graph, &partitions, &clock_model)?;

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
      .ok_or_else(|| eyre::eyre!("--sequence-length required for --covariation without alignment"))? as f64
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

fn write_outputs(
  args: &TreetimeTimetreeArgs,
  graph: &GraphTimetree,
  _partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  clock_model: &ClockModel,
) -> Result<(), Report> {
  let out_base = args.outdir.join("timetree");
  nwk_write_file(out_base.with_extension("nwk"), graph, &NwkWriteOptions::default())
    .wrap_err("Failed to write Newick output")?;
  nex_write_file(out_base.with_extension("nexus"), graph, &NexWriteOptions::default())
    .wrap_err("Failed to write Nexus output")?;

  write_clock_model(clock_model, &out_base)?;

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
