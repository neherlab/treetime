use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::convergence::metrics::{IterationContext, TimetreeOptimizer};
use crate::commands::timetree::inference::runner::run_timetree;
use crate::commands::timetree::initialization::{InputData, load_input_data};
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::commands::timetree::refinement::run_refinement_iteration;
use crate::representation::graph_ancestral::GraphAncestral;
use eyre::{Report, WrapErr};
use log::info;
use parking_lot::RwLock;
use std::sync::Arc;

pub fn run_timetree_estimation(args: &TreetimeTimetreeArgs) -> Result<(), Report> {
  log::info!("# TreeTime Timetree Estimation");
  log::debug!("Branch length mode: {:?}, Keep root: {}", args.branch_length_mode, args.keep_root);

  log::info!("## Loading input data");
  let InputData { mut graph, aln, .. } = load_input_data(args)?;
  log::info!("Input loaded: {} nodes, {} edges", graph.get_nodes().len(), graph.get_edges().len());

  // let mut clock_model = infer_clock_model(args, &graph).wrap_err("Failed to infer clock model")?;

  // In input branch mode, skip partition initialization - we use branch lengths from the tree directly
  let partitions = vec![];

  // run_pre_optimization(args, &graph, &partitions, aln.as_deref())?;

  run_initial_timetree_inference(args, &mut graph, &partitions)?;

  let mut optimizer = TimetreeOptimizer::new(args.max_iter, args.tracelog.clone())?;
  while let Some(IterationContext { i }) = optimizer.next_iter() {
    let (ndiff, n_resolved) = run_refinement_iteration(args, &mut graph, &partitions, aln.as_deref(), i)?;

    // clock_model = update_clock_model(&graph, &clock_model)
    // .wrap_err_with(|| format!("Failed to update clock model (iteration {i})"))?;

    optimizer
      .record(ndiff, n_resolved, &partitions)
      .wrap_err_with(|| format!("Failed to record convergence metrics (iteration {i})"))?;
  }

  // run_post_processing(args, &graph, &partitions, &clock_model)?;

  // write_outputs(args, &graph, &partitions, &clock_model)?;

  Ok(())
}

// fn run_pre_optimization(
//   args: &TreetimeTimetreeArgs,
//   graph: &GraphAncestral,
//   partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
//   aln: Option<&[FastaRecord]>,
// ) -> Result<(), Report> {
//   // if !args.keep_root {
//   //   info!("Rerooting tree to optimize temporal signal");
//   //   todo!("reroot_tree not yet updated for GraphAncestral");
//   // }

//   // if args.clock_filter_enabled() {
//   //   todo!("clock_filter not yet updated for GraphAncestral");
//   // }

//   if aln.is_some() {
//     match args.branch_length_mode {
//       BranchLengthMode::Input => {
//         info!("Using input branch lengths for timetree inference");
//       },
//       BranchLengthMode::Marginal | BranchLengthMode::Auto => {
//         info!("Running initial ancestral reconstruction");
//         // For sparse: compress_sequences populates partitions via Fitch, then run_marginal does message passing
//         // For dense: run_marginal with aln argument handles initialization via attach_sequences
//         // Since partitions are trait objects, we can't call compress_sequences directly
//         // Instead, run_marginal will handle both cases appropriately
//         run_marginal(graph, partitions, aln)?;
//       },
//       BranchLengthMode::Joint => {
//         todo!("Joint ancestral reconstruction not yet implemented")
//       },
//     }
//   }

//   Ok(())
// }

fn run_initial_timetree_inference(
  args: &TreetimeTimetreeArgs,
  graph: &mut GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  // _aln: Option<&[FastaRecord]>,
  // _clock_model: &ClockModel,
) -> Result<(), Report> {
  info!("### TreeTime: INITIAL ROUND");
  run_timetree(graph, partitions, args.keep_root)?;

  // if !_args.keep_root {
  //   todo!("Second reroot not yet implemented");
  // }

  Ok(())
}

// fn run_post_processing(
//   args: &TreetimeTimetreeArgs,
//   graph: &GraphAncestral,
//   partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
//   _clock_model: &ClockModel,
// ) -> Result<(), Report> {
//   if args.vary_rate {
//     todo!("calc_rate_susceptibility not yet implemented");
//   }

//   if args.time_marginal == TimeMarginalMode::OnlyFinal {
//     info!("### Final round: marginal reconstruction for confidence intervals");
//     run_timetree(graph, partitions).wrap_err("Final timetree inference failed")?;
//     todo!("extract_confidence_intervals not yet implemented");
//   }

//   Ok(())
// }

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
