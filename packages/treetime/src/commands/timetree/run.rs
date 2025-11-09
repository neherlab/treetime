use crate::commands::ancestral::marginal_unified::run_marginal;
use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_regression::{ClockOptions, estimate_clock_model_with_reroot};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::timetree::args::{BranchLengthMode, TimeMarginalMode, TreetimeTimetreeArgs};
use crate::commands::timetree::convergence::metrics::{IterationContext, TimetreeOptimizer};
use crate::commands::timetree::inference::runner::run_timetree;
use crate::commands::timetree::initialization::{InputData, initialize_partitions, load_input_data};
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::commands::timetree::refinement::run_refinement_iteration;
use crate::commands::timetree::utils::initialize_clock_totals_from_time_distributions;
use crate::io::fasta::FastaRecord;
use crate::representation::graph_ancestral::GraphAncestral;
use eyre::{Report, WrapErr};
use log::info;
use parking_lot::RwLock;
use std::sync::Arc;

pub fn run_timetree_estimation(args: &TreetimeTimetreeArgs) -> Result<(), Report> {
  log::info!("# TreeTime Timetree Estimation");
  log::debug!(
    "Branch length mode: {:?}, Keep root: {}",
    args.branch_length_mode,
    args.keep_root
  );

  log::info!("## Loading input data");
  let InputData {
    mut graph,
    alphabet,
    aln,
  } = load_input_data(args)?;
  log::info!(
    "Input loaded: {} nodes, {} edges",
    graph.get_nodes().len(),
    graph.get_edges().len()
  );

  let clock_options = ClockOptions::default();
  let branch_params = BranchPointOptimizationParams::default();

  let mut clock_model = estimate_clock_model_with_reroot(
    &mut graph,
    &clock_options,
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
      initialize_partitions(args, &graph, alphabet, aln.as_deref())?
    },
  };

  run_pre_optimization(args, &graph, &partitions, aln.as_deref())?;

  run_initial_timetree_inference(args, &mut graph, &partitions, &clock_model)?;

  let mut optimizer = TimetreeOptimizer::new(args.max_iter, args.tracelog.clone())?;
  while let Some(IterationContext { i }) = optimizer.next_iter() {
    let (ndiff, n_resolved) = run_refinement_iteration(args, &mut graph, &partitions, aln.as_deref(), i, &clock_model)?;

    clock_model = estimate_clock_model_with_reroot(
      &mut graph,
      &clock_options,
      args.clock_rate,
      args.keep_root,
      &branch_params,
    )
    .wrap_err_with(|| format!("Failed to update clock model (iteration {i})"))?;

    optimizer
      .record(ndiff, n_resolved, &partitions)
      .wrap_err_with(|| format!("Failed to record convergence metrics (iteration {i})"))?;
  }

  run_post_processing(args, &mut graph, &partitions, &clock_model)?;

  // write_outputs(args, &graph, &partitions, &clock_model)?;

  Ok(())
}

fn run_pre_optimization(
  args: &TreetimeTimetreeArgs,
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  aln: Option<&[FastaRecord]>,
) -> Result<(), Report> {
  // if !args.keep_root {
  //   info!("Rerooting tree to optimize temporal signal");
  //   todo!("reroot_tree not yet updated for GraphAncestral");
  // }

  // if args.clock_filter_enabled() {
  //   todo!("clock_filter not yet updated for GraphAncestral");
  // }

  if aln.is_some() {
    match args.branch_length_mode {
      BranchLengthMode::Input => {
        info!("Using input branch lengths for timetree inference");
      },
      BranchLengthMode::Marginal => {
        info!("Running initial ancestral reconstruction");
        run_marginal(graph, partitions, aln)?;
      },
    }
  }

  Ok(())
}

fn run_initial_timetree_inference(
  args: &TreetimeTimetreeArgs,
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  clock_model: &ClockModel,
) -> Result<(), Report> {
  info!("### TreeTime: INITIAL ROUND");
  info!("### Initializing node times from date constraints");
  initialize_clock_totals_from_time_distributions(graph)?;

  run_timetree(graph, partitions, clock_model)?;

  // if !_args.keep_root {
  //   todo!("Second reroot not yet implemented");
  // }

  Ok(())
}

fn run_post_processing(
  args: &TreetimeTimetreeArgs,
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  clock_model: &ClockModel,
) -> Result<(), Report> {
  if args.vary_rate {
    todo!("calc_rate_susceptibility not yet implemented");
  }

  if args.time_marginal == TimeMarginalMode::OnlyFinal {
    info!("### Final round: marginal reconstruction for confidence intervals");
    run_timetree(graph, partitions, clock_model).wrap_err("Final timetree inference failed")?;
    todo!("extract_confidence_intervals not yet implemented");
  }

  Ok(())
}

// fn write_outputs(
//   args: &TreetimeTimetreeArgs,
//   graph: &GraphAncestral,
//   _partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
//   clock_model: &ClockModel,
// ) -> Result<(), Report> {
//   let out_base = args.outdir.join("timetree");
//   nwk_write_file(out_base.with_extension("nwk"), graph, &NwkWriteOptions::default())
//     .wrap_err("Failed to write Newick output")?;
//   nex_write_file(out_base.with_extension("nexus"), graph, &NexWriteOptions::default())
//     .wrap_err("Failed to write Nexus output")?;

//   write_clock_model(clock_model, &out_base)?;

//   if args.plot_rtt.is_some() {
//     todo!("plot_root_to_tip not yet updated for GraphAncestral");
//   }

//   if args.plot_tree.is_some() {
//     todo!("plot_time_tree not yet updated for GraphAncestral");
//   }

//   Ok(())
// }
