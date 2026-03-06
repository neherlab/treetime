use crate::commands::ancestral::marginal::initialize_marginal;
use crate::commands::clock::clock_filter::clock_filter_inplace;
use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_output::write_clock_model;
use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::timetree::args::{BranchLengthMode, TimeMarginalMode, TreetimeTimetreeArgs};
use crate::commands::timetree::coalescent::optimize_tc::optimize_tc;
use crate::commands::timetree::coalescent::skyline::{SkylineParams, optimize_skyline};
use crate::commands::timetree::convergence::metrics::{IterationContext, TimetreeOptimizer};
use crate::commands::timetree::inference::runner::run_timetree;
use crate::commands::timetree::initialization::{InputData, initialize_partitions, load_input_data};
use crate::commands::timetree::optimization::clock_filter::report_bad_branches;
use crate::commands::timetree::optimization::reroot::reroot_tree;
use crate::commands::timetree::output::confidence::{extract_confidence_intervals, write_confidence_intervals};
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::commands::timetree::refinement::run_refinement_iteration;
use crate::commands::timetree::utils::initialize_clock_totals_from_time_distributions;
use crate::gtr::get_gtr::{GtrModelName, JC69Params, jc69, write_gtr_json};
use crate::representation::partition::timetree::GraphTimetree;
use crate::representation::payload::timetree::EdgeTimetree;
use crate::representation::payload::timetree::NodeTimetree;
use eyre::{Report, WrapErr};
use log::{debug, info, warn};
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_io::nex::{NexWriteOptions, nex_write_file};
use treetime_io::nwk::{NwkWriteOptions, nwk_write_file};

pub fn run_timetree_estimation(args: &TreetimeTimetreeArgs) -> Result<(), Report> {
  info!("# TreeTime Timetree Estimation");
  debug!(
    "Branch length mode: {:?}, Keep root: {}",
    args.branch_length_mode, args.keep_root
  );

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

  let clock_params = ClockParams::default();
  let branch_params = BranchPointOptimizationParams::default();

  let mut clock_model = estimate_clock_model_with_reroot(
    &mut graph,
    &clock_params,
    args.clock_rate,
    args.keep_root,
    &branch_params,
  )
  .wrap_err("Failed to infer clock model")?;

  let partitions = match args.branch_length_mode {
    BranchLengthMode::Input => {
      info!("Branch length mode: Input - using tree branch lengths");
      vec![]
    },
    BranchLengthMode::Marginal => {
      info!("Branch length mode: Marginal - initializing partitions from alignment");
      let partitions = initialize_partitions(args, &graph, alphabet, aln.as_deref())?;
      // Write GTR model used for marginal reconstruction
      // FIXME: currently hardcoded to JC69, should use args.gtr when model selection is implemented
      let gtr = jc69(JC69Params::default())?;
      write_gtr_json(&gtr, GtrModelName::JC69, &args.outdir)?;
      partitions
    },
  };

  if !args.keep_root {
    info!("First reroot (pre-ancestral)");
    clock_model = reroot_tree(&mut graph, &partitions, &clock_params, args.clock_rate, &branch_params)
      .wrap_err("Failed to reroot tree (pre-ancestral)")?;
  }

  if args.clock_filter > 0.0 {
    let result = clock_filter_inplace(&graph, &clock_model, args.clock_filter);
    report_bad_branches(&graph, &clock_model, result.iqd);
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

  // Create coalescent Tc distribution
  let mut coalescent_tc: Option<Distribution> = if args.coalescent_skyline {
    info!(
      "### Optimizing skyline coalescent model with {} grid points",
      args.n_skyline
    );
    let skyline_result =
      optimize_skyline(&graph, &skyline_params).wrap_err("Failed to optimize skyline coalescent model")?;
    info!(
      "Skyline optimization completed: log_likelihood={:.4}",
      skyline_result.log_likelihood
    );
    Some(skyline_result.tc_distribution)
  } else {
    args.coalescent.map(Distribution::constant)
  };

  run_timetree(&mut graph, &partitions, &clock_model, coalescent_tc.as_ref())?;

  if !args.keep_root {
    info!("Reroot (post-ancestral)");
    clock_model = reroot_tree(&mut graph, &partitions, &clock_params, args.clock_rate, &branch_params)
      .wrap_err("Failed to reroot tree (post-ancestral)")?;
  }

  info!("### TreeTime: Optimisation rounds");
  let mut optimizer = TimetreeOptimizer::new(args.max_iter, args.tracelog.clone())?;
  while let Some(IterationContext { i }) = optimizer.next_iter() {
    // Optimize Tc if requested (requires at least 2 iterations per Python v0)
    // Note: Tc optimization only applies when not using skyline (skyline has its own optimization)
    if args.coalescent_opt && !args.coalescent_skyline && i >= 2 {
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
      &clock_params,
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

    // Run final timetree pass with optimized skyline, unless OnlyFinal mode will run it below
    if args.time_marginal != TimeMarginalMode::OnlyFinal {
      run_timetree(&mut graph, &partitions, &clock_model, coalescent_tc.as_ref())
        .wrap_err("Final timetree pass with optimized skyline failed")?;
    }
  }

  info!("### TreeTime: postprocessing");
  if args.vary_rate {
    todo!("calc_rate_susceptibility not yet implemented");
  }

  if args.time_marginal == TimeMarginalMode::OnlyFinal {
    info!("### Final round: marginal reconstruction for confidence intervals");
    run_timetree(&mut graph, &partitions, &clock_model, coalescent_tc.as_ref())
      .wrap_err("Final timetree inference failed")?;
    let intervals = extract_confidence_intervals(&graph);
    let ci_path = args.outdir.join("confidence_intervals.tsv");
    write_confidence_intervals(&intervals, &ci_path).wrap_err("Failed to write confidence intervals")?;
    info!("Wrote confidence intervals to {ci_path}", ci_path = ci_path.display());
  }

  info!("### TreeTime: writing outputs");
  write_outputs(args, &graph, &partitions, &clock_model)?;

  Ok(())
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
    todo!("plot_root_to_tip not yet implemented");
  }

  if args.plot_tree.is_some() {
    todo!("plot_time_tree not yet implemented");
  }

  Ok(())
}
