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
use crate::partition::timetree::{GraphTimetree, PartitionTimetreeAllVec};
use crate::partition::traits::{HasGtr, PartitionTimetreeAll};
use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
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
const INITIAL_COALESCENT_TC: f64 = 5.0;

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
    ..SkylineParams::default()
  };

  if params.n_branches_posterior.is_some() {
    return make_error!("--n-branches-posterior is not yet implemented");
  }

  run_timetree(&mut input.graph, &partitions, &clock_model, None, params.no_indels)?;

  let mut coalescent_tc: Option<Distribution> =
    match coalescent_initialization(params.coalescent, params.coalescent_opt, params.coalescent_skyline) {
      CoalescentInitialization::Optimize(initial_tc) => {
        info!("### Optimizing constant coalescent Tc before refinement");
        match optimize_tc(&input.graph, initial_tc) {
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
      },
      CoalescentInitialization::Fixed(tc) => Some(Distribution::constant(tc)),
      CoalescentInitialization::None => None,
    };

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
  let mut optimizer = TimetreeOptimizer::new(params.max_iter, params.coalescent_skyline);
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

    if (params.coalescent_opt || params.coalescent_skyline) && i >= 2 {
      let initial_tc = coalescent_tc.as_ref().map_or(1.0, |d| d.max_value());
      match optimize_tc(&input.graph, initial_tc) {
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

  if params.coalescent_skyline {
    info!("### Re-optimizing skyline coalescent with stabilized node times");
    let skyline_result =
      optimize_skyline(&input.graph, &skyline_params).wrap_err("Failed to re-optimize skyline coalescent model")?;
    info!(
      "Skyline re-optimization completed: log_likelihood={:.4}",
      skyline_result.log_likelihood
    );
    coalescent_tc = Some(skyline_result.tc_distribution);

    if time_marginal != TimeMarginalMode::OnlyFinal {
      run_timetree(
        &mut input.graph,
        &partitions,
        &clock_model,
        coalescent_tc.as_ref(),
        params.no_indels,
      )
      .wrap_err("Final timetree pass with optimized skyline failed")?;

      if !partitions.is_empty() {
        update_marginal(&input.graph, &partitions)?;
      }
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

#[derive(Clone, Copy, Debug, PartialEq)]
enum CoalescentInitialization {
  Optimize(f64),
  Fixed(f64),
  None,
}

fn coalescent_initialization(
  coalescent: Option<f64>,
  coalescent_opt: bool,
  coalescent_skyline: bool,
) -> CoalescentInitialization {
  if coalescent_opt || coalescent_skyline {
    CoalescentInitialization::Optimize(coalescent.unwrap_or(INITIAL_COALESCENT_TC))
  } else if let Some(tc) = coalescent {
    CoalescentInitialization::Fixed(tc)
  } else {
    CoalescentInitialization::None
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

  let partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> = match created.partition {
    MarginalPartition::Sparse(p) => Arc::new(RwLock::new(p)),
    MarginalPartition::Dense(p) => Arc::new(RwLock::new(p)),
  };

  Ok(PartitionInitResult {
    partitions: vec![partition],
    gtr,
    model_name: created.model_name,
  })
}

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

  apply_damping(graph, &old_branch_lengths, TIMETREE_PRE_STEP_DAMPING, 0);
  update_marginal(graph, partitions)?;

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::{CoalescentInitialization, coalescent_initialization};
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::disabled(       None,       false, false, CoalescentInitialization::None)]
  #[case::fixed(          Some(0.25), false, false, CoalescentInitialization::Fixed(0.25))]
  #[case::opt_default(    None,       true,  false, CoalescentInitialization::Optimize(0.001))]
  #[case::opt_configured( Some(0.25), true,  false, CoalescentInitialization::Optimize(0.25))]
  #[case::skyline_default(None,       false, true,  CoalescentInitialization::Optimize(0.001))]
  #[trace]
  fn test_pipeline_coalescent_initialization(
    #[case] coalescent: Option<f64>,
    #[case] coalescent_opt: bool,
    #[case] coalescent_skyline: bool,
    #[case] expected: CoalescentInitialization,
  ) {
    let actual = coalescent_initialization(coalescent, coalescent_opt, coalescent_skyline);

    // v0 creates Coalescent with Tc=0.001 before optimizing it in the first iteration:
    // packages/legacy/treetime/treetime/merger_models.py#L25
    // packages/legacy/treetime/treetime/treetime.py#L307-L317
    assert_eq!(expected, actual);
  }
}
