use crate::commands::ancestral::marginal_unified::run_marginal;
use crate::commands::timetree::args::{BranchLengthMode, TimeMarginalMode, TreetimeTimetreeArgs};
use crate::commands::timetree::convergence::metrics::{IterationContext, TimetreeOptimizer};
use crate::commands::timetree::data::clock_model::{ClockModel, infer_clock_model, update_clock_model};
use crate::commands::timetree::data::date_constraints::DateConstraintSet;
use crate::commands::timetree::inference::runner::run_timetree;
use crate::commands::timetree::initialization::{InputData, initialize_partitions, load_input_data};
use crate::commands::timetree::optimization::clock_filter::{clock_filter, report_bad_branches};
use crate::commands::timetree::optimization::reroot::reroot_tree;
use crate::commands::timetree::output::clock_model_output::write_clock_model;
use crate::commands::timetree::output::confidence::{
  calc_rate_susceptibility, extract_confidence_intervals, write_confidence_intervals,
};
use crate::commands::timetree::output::dates::write_node_dates;
use crate::commands::timetree::output::plots::{plot_root_to_tip, plot_time_tree};
use crate::commands::timetree::refinement::run_refinement_iteration;
use crate::io::fasta::FastaRecord;
use crate::io::nex::{NexWriteOptions, nex_write_file};
use crate::io::nwk::{NwkWriteOptions, nwk_write_file};
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::partition_timetree::PartitionTreetimeMarginalOps;
use eyre::{Report, WrapErr};
use log::info;
use parking_lot::RwLock;
use std::sync::Arc;

pub fn run_timetree_estimation(args: &TreetimeTimetreeArgs) -> Result<(), Report> {
  let InputData {
    graph,
    alphabet,
    aln,
    constraints,
  } = load_input_data(args)?;

  let mut clock_model = infer_clock_model(args, &graph, &constraints).wrap_err("Failed to infer clock model")?;

  let partitions: Vec<Arc<RwLock<dyn PartitionTreetimeMarginalOps>>> =
    initialize_partitions(args, &graph, alphabet, aln.as_deref(), &constraints)?;

  run_pre_optimization(args, &graph, &partitions, aln.as_deref(), &constraints)?;

  run_initial_timetree_inference(args, &graph, &partitions, aln.as_deref(), &constraints)?;

  let mut optimizer = TimetreeOptimizer::new(args.max_iter, args.tracelog.clone())?;
  while let Some(IterationContext { i }) = optimizer.next_iter() {
    let (ndiff, n_resolved) = run_refinement_iteration(args, &graph, &partitions, aln.as_deref(), i)?;

    clock_model = update_clock_model(&graph, &constraints, &clock_model)
      .wrap_err_with(|| format!("Failed to update clock model (iteration {i})"))?;

    optimizer
      .record(ndiff, n_resolved, &partitions)
      .wrap_err_with(|| format!("Failed to record convergence metrics (iteration {i})"))?;
  }

  run_post_processing(args, &graph, &partitions, &constraints, &clock_model)?;

  write_outputs(args, &graph, &partitions, &constraints, &clock_model)?;

  Ok(())
}

fn run_pre_optimization(
  args: &TreetimeTimetreeArgs,
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
  aln: Option<&[FastaRecord]>,
  constraints: &DateConstraintSet,
) -> Result<(), Report> {
  if !args.keep_root {
    info!("Rerooting tree to optimize temporal signal");
    reroot_tree(graph, constraints, "least-squares")?;
  }

  // TODO: Implement clock filter functionality
  // What: Detect and exclude branches that violate molecular clock assumptions using IQD threshold
  // Why: Optional - improves robustness by removing outliers from clock model inference
  // How: Analyze root-to-tip regression residuals; mark branches exceeding n_iqd * IQD as bad_branch=true
  if args.clock_filter_enabled() {
    clock_filter(graph, constraints, args.n_iqd.unwrap_or(3.0))?;
  }

  // Initial ancestral sequence reconstruction (establishes branch length distributions)
  if aln.is_some() {
    match args.branch_length_mode {
      BranchLengthMode::Input => {
        info!("Using input branch lengths for timetree inference");
      },
      BranchLengthMode::Marginal | BranchLengthMode::Auto => {
        info!("Running initial ancestral reconstruction");
        run_marginal(graph, partitions, aln)?;
      },
      BranchLengthMode::Joint => {
        todo!("Joint ancestral reconstruction not yet implemented")
      },
    }
  }

  Ok(())
}

fn run_initial_timetree_inference(
  args: &TreetimeTimetreeArgs,
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
  aln: Option<&[FastaRecord]>,
  constraints: &DateConstraintSet,
) -> Result<(), Report> {
  info!("### TreeTime: INITIAL ROUND");
  run_timetree(graph, partitions)?;

  // TODO: Implement second rerooting with inferred node times
  // What: Reroot tree again using node times inferred from initial timetree run
  // Why: Optional - initial rerooting used sampling dates only; second pass uses better time estimates
  // How: Run reroot with inferred node times, then re-run ancestral reconstruction and timetree inference
  if !args.keep_root {
    reroot_tree(graph, constraints, "least-squares")?;
    if args.branch_length_mode != BranchLengthMode::Input {
      run_marginal(graph, partitions, aln)?;
    }
    run_timetree(graph, partitions)?;
  }

  Ok(())
}

fn run_post_processing(
  args: &TreetimeTimetreeArgs,
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
  constraints: &DateConstraintSet,
  clock_model: &ClockModel,
) -> Result<(), Report> {
  // TODO: Implement rate variation sensitivity analysis
  // What: Re-run timetree inference with rate ± standard deviation to quantify uncertainty propagation
  // Why: Optional - assesses robustness of time estimates to clock rate uncertainty
  // How: Run timetree with perturbed rates, store alternative node time estimates for comparison
  if args.vary_rate {
    calc_rate_susceptibility(graph, partitions, constraints, clock_model)?;
  }

  // Final marginal reconstruction for confidence intervals
  if args.time_marginal == TimeMarginalMode::OnlyFinal {
    info!("### Final round: marginal reconstruction for confidence intervals");
    run_timetree(graph, partitions).wrap_err("Final timetree inference failed")?;
    // TODO: Implement confidence interval extraction
    // What: Extract 95% HPD intervals from marginal posterior distributions for each node
    // Why: Core - quantifies uncertainty in inferred divergence times
    // How: Compute highest posterior density regions or quantiles from marginal_pos_LH distributions
    extract_confidence_intervals(graph, partitions)?;
  }

  // Identify and report outlier branches
  // TODO: Implement bad branch reporting
  // What: Identify branches where inferred dates deviate significantly from input constraints
  // Why: Core - users need to know which samples were excluded or flagged as problematic
  // How: Compare inferred dates to input constraints, log warnings for large deviations
  report_bad_branches(graph, constraints)?;

  Ok(())
}

fn write_outputs(
  args: &TreetimeTimetreeArgs,
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
  constraints: &DateConstraintSet,
  clock_model: &ClockModel,
) -> Result<(), Report> {
  let out_base = args.outdir.join("timetree");
  nwk_write_file(out_base.with_extension("nwk"), graph, &NwkWriteOptions::default())
    .wrap_err("Failed to write Newick output")?;
  nex_write_file(out_base.with_extension("nexus"), graph, &NexWriteOptions::default())
    .wrap_err("Failed to write Nexus output")?;

  // TODO: Implement node dates output
  // What: Write TSV file with inferred calendar dates for each tree node
  // Why: Core - primary output showing numdate and time_before_present for each node
  // How: Extract node time from partition data, convert to calendar dates, write TSV
  write_node_dates(graph, partitions, &out_base)?;

  // TODO: Implement confidence intervals output
  // What: Write TSV file with lower_bound, median, upper_bound for each node's inferred date
  // Why: Core - quantifies uncertainty in molecular dating results
  // How: Extract HPD intervals from marginal distributions, write TSV with node_name and bounds
  write_confidence_intervals(graph, partitions, &out_base)?;

  // TODO: Implement clock model output
  // What: Write JSON file with clock rate, intercept, R², and confidence intervals
  // Why: Core - documents molecular clock parameters used in analysis
  // How: Serialize clock_model struct to JSON with rate, std_dev, intercept, R² fields
  write_clock_model(clock_model, &out_base)?;

  // TODO: Implement root-to-tip regression plot
  // What: Generate scatter plot of root-to-tip distance vs sampling date with regression line
  // Why: Optional - diagnostic visualization for temporal signal quality and outlier detection
  // How: Extract leaf dates and distances, plot with matplotlib/plotters, annotate with R² and outliers
  if args.plot_rtt.is_some() {
    plot_root_to_tip(graph, constraints, &out_base)?;
  }

  // TODO: Implement time-scaled tree plot
  // What: Generate visualization of phylogenetic tree with x-axis as calendar time
  // Why: Optional - helps interpret divergence times and tree topology in temporal context
  // How: Layout tree horizontally with branch lengths scaled by inferred divergence times
  if args.plot_tree.is_some() {
    plot_time_tree(graph, &out_base)?;
  }

  Ok(())
}
