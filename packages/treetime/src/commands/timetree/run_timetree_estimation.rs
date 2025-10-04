#[allow(dead_code, clippy::wildcard_imports)]
use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::get_common_length;
use crate::commands::ancestral::marginal_unified::run_marginal;
use crate::commands::timetree::clock_model::ClockModel;
use crate::commands::timetree::clock_model::{infer_clock_model, update_clock_model};
use crate::commands::timetree::date_constraints::DateConstraintSet;
use crate::commands::timetree::date_constraints::load_date_constraints;
use crate::commands::timetree::timetree_args::{BranchLengthMode, TimeMarginalMode, TreetimeTimetreeArgs};
use crate::commands::timetree::timetree_unified::run_timetree;
use crate::gtr::get_gtr::{JC69Params, jc69};
use crate::io::csv::CsvStructFileWriter;
use crate::io::fasta::{FastaRecord, read_many_fasta};
use crate::io::nex::{NexWriteOptions, nex_write_file};
use crate::io::nwk::{NwkWriteOptions, nwk_read_file, nwk_write_file};
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::infer_dense::infer_dense;
use crate::representation::partition_timetree::PartitionTimetreeOps;
use crate::representation::partition_timetree::PartitionTreetimeMarginalOps;
use crate::representation::partition_timetree_dense::PartitionTimetreeDense;
use crate::representation::partition_timetree_sparse::PartitionTimetreeSparse;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::info;
use maplit::btreemap;
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

// ============================================================================
// Convergence tracking data structures
// ============================================================================

/// Tracks convergence metrics across timetree optimization iterations.
///
/// Records likelihood components and change counts to monitor convergence:
/// - Sequence changes (ndiff) should approach zero
/// - Polytomies resolved (n_resolved) should stabilize
/// - Likelihoods should increase or stabilize
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConvergenceMetrics {
  /// Number of ancestral sequence changes in this iteration
  pub ndiff: usize,
  /// Number of polytomies resolved in this iteration
  pub n_resolved: usize,
  /// Sequence likelihood (probability of observing sequences given tree and substitution model)
  pub lh_seq: Option<f64>,
  /// Positional likelihood (probability of node positions on time axis)
  pub lh_pos: Option<f64>,
  /// Coalescent likelihood (population genetic prior on node times)
  pub lh_coal: Option<f64>,
  /// Total likelihood (combined probability of all components)
  pub lh_total: Option<f64>,
}

impl ConvergenceMetrics {
  pub fn from_iteration(
    ndiff: usize,
    n_resolved: usize,
    partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
  ) -> Self {
    Self {
      ndiff,
      n_resolved,
      lh_seq: calculate_sequence_likelihood(partitions),
      lh_pos: calculate_positional_likelihood(partitions),
      lh_coal: calculate_coalescent_likelihood(partitions),
      lh_total: calculate_total_likelihood(partitions),
    }
  }

  pub fn has_converged(&self) -> bool {
    self.ndiff == 0 && self.n_resolved == 0
  }
}

struct TreetimeOptimizerTraceCsvWriter {
  writer: CsvStructFileWriter,
}

impl TreetimeOptimizerTraceCsvWriter {
  fn new(path: impl AsRef<Path>) -> Result<Self, Report> {
    let writer = CsvStructFileWriter::new(path, b',')?;
    Ok(Self { writer })
  }

  fn write(&mut self, metrics: &ConvergenceMetrics) -> Result<(), Report> {
    self.writer.write(metrics)?;
    Ok(())
  }
}

struct TimetreeOptimizer {
  trace: Vec<ConvergenceMetrics>,
  tracelog_writer: Option<TreetimeOptimizerTraceCsvWriter>,
  max_iterations: usize,
  i: usize,
}

impl TimetreeOptimizer {
  fn new(max_iter: usize, tracelog_path: Option<impl Into<PathBuf>>) -> Result<Self, Report> {
    let tracelog_writer = tracelog_path
      .map(Into::into)
      .map(TreetimeOptimizerTraceCsvWriter::new)
      .transpose()?;

    Ok(Self {
      trace: vec![],
      tracelog_writer,
      max_iterations: max_iter,
      i: 0,
    })
  }

  fn record(
    &mut self,
    n_diff: usize,
    n_resolved: usize,
    partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
  ) -> Result<(), Report> {
    let metric = ConvergenceMetrics::from_iteration(n_diff, n_resolved, partitions);

    if let Some(writer) = &mut self.tracelog_writer {
      writer.write(&metric)?;
    }

    if metric.has_converged() {
      info!(
        "Converged at iteration {} (ndiff={n_diff}, n_resolved={n_resolved})",
        self.i
      );
    } else {
      info!(
        "  Iteration {}: n_diff={n_diff}, n_resolved={n_resolved}, total_LH={:.2}",
        self.i,
        metric.lh_total.unwrap_or(f64::NAN)
      );
    }

    self.trace.push(metric);
    Ok(())
  }

  fn has_converged(&self) -> bool {
    self.trace.last().is_some_and(|m| m.has_converged())
  }

  fn has_reached_max_iterations(&self) -> bool {
    self.i >= self.max_iterations
  }

  fn next_iter(&mut self) -> Option<IterationContext> {
    if self.has_converged() || self.has_reached_max_iterations() {
      return None;
    }

    self.i += 1;
    info!("### Timetree iteration {}/{}", self.i, self.max_iterations);

    Some(IterationContext { i: self.i })
  }
}

pub struct IterationContext {
  i: usize,
}

/// Stores ancestral sequence states for change tracking between iterations.
///
/// Maps node keys to their reconstructed sequences to enable ndiff calculation.
type AncestralStateSnapshot = BTreeMap<String, Vec<u8>>;

struct InputData {
  graph: GraphAncestral,
  alphabet: Alphabet,
  aln: Option<Vec<FastaRecord>>,
  constraints: DateConstraintSet,
}

fn load_input_data(args: &TreetimeTimetreeArgs) -> Result<InputData, Report> {
  let graph: GraphAncestral = if let Some(tree_path) = &args.tree {
    nwk_read_file(tree_path).wrap_err("Failed to load tree from file")?
  } else {
    todo!("Tree inference from alignment not yet implemented")
  };

  // Create alphabet for both alignment loading and timetree inference.
  // treat_gap_as_unknown behavior:
  // - With alignment: uses dense setting (true for dense mode, false for sparse)
  // - Without alignment: uses false (gaps as distinct characters for branch length mode)
  let alphabet = {
    let treat_gap_as_unknown = !args.input_fastas.is_empty() && args.dense.unwrap_or_else(infer_dense);
    Alphabet::new(args.alphabet, treat_gap_as_unknown)?
  };

  // Load alignment sequences (optional if using input branch lengths only)
  let aln = if !args.input_fastas.is_empty() {
    Some(read_many_fasta(&args.input_fastas, &alphabet)?)
  } else if args.input_fastas.is_empty() && args.branch_length_mode != BranchLengthMode::Input {
    return Err(eyre::eyre!(
      "Alignment required when branch_length_mode is not 'input'. \
       Provide FASTA files or use --branch-length-mode=input"
    ));
  } else {
    None
  };

  let constraints = load_date_constraints(args, &graph).wrap_err("Failed to load date constraints")?;

  Ok(InputData {
    graph,
    alphabet,
    aln,
    constraints,
  })
}

fn initialize_partitions(
  args: &TreetimeTimetreeArgs,
  graph: &GraphAncestral,
  alphabet: Alphabet,
  aln: Option<&[FastaRecord]>,
  constraints: &DateConstraintSet,
) -> Result<Vec<Arc<RwLock<dyn PartitionTreetimeMarginalOps>>>, Report> {
  let dense = args.dense.unwrap_or_else(infer_dense);
  let sequence_length = if let Some(aln_data) = aln {
    Some(get_common_length(aln_data)?)
  } else {
    args.sequence_length
  };

  #[allow(clippy::iter_on_single_items, trivial_casts)]
  let partitions = if !dense {
    [Arc::new(RwLock::new(PartitionTimetreeSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      sequence_length,
      nodes: btreemap! {},
      edges: btreemap! {},
    })) as Arc<RwLock<dyn PartitionTreetimeMarginalOps>>]
  } else {
    [Arc::new(RwLock::new(PartitionTimetreeDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      sequence_length,
      nodes: btreemap! {},
      edges: btreemap! {},
    })) as Arc<RwLock<dyn PartitionTreetimeMarginalOps>>]
  }
  .into_iter()
  .map(|partition| -> Result<_, Report> {
    PartitionTimetreeOps::attach_date_constraints(&mut *partition.write_arc(), graph, &constraints.per_node)?;
    Ok(partition)
  })
  .try_collect()?;

  Ok(partitions)
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

#[allow(clippy::useless_let_if_seq)]
fn run_refinement_iteration(
  args: &TreetimeTimetreeArgs,
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
  aln: Option<&[FastaRecord]>,
  i: usize,
) -> Result<(usize, usize), Report> {
  let mut is_tree_dirty = false;

  // Add coalescent model (population genetic prior on node times)
  // TODO: Implement coalescent model integration
  // What: Add Kingman coalescent prior to node time distributions using population genetic theory
  // Why: Optional - incorporates realistic population dynamics to improve time estimates
  // How: Calculate merger rates between lineages, add to node time likelihood as prior term
  if let Some(coalescent_params) = &args.coalescent {
    add_coalescent_model(graph, partitions, coalescent_params)?;
    is_tree_dirty = true;
  }

  // Apply relaxed clock (allow branch-specific rate variation)
  // TODO: Implement relaxed clock model
  // What: Allow each branch to have its own rate multiplier (gamma) with autocorrelation constraints
  // Why: Optional - accounts for rate heterogeneity when strict clock assumption is violated
  // How: Optimize branch-specific gamma values with slack/coupling penalties to prevent overfitting
  if !args.relax.is_empty() {
    apply_relaxed_clock(graph, partitions, &args.relax)?;
    is_tree_dirty = true;
  }

  // Resolve polytomies using temporal constraints
  // TODO: Implement polytomy resolution algorithm
  // What: Convert multifurcations into binary trees using likelihood-based greedy optimization
  // Why: Optional - improves tree topology and is required for some downstream analyses
  // How: For each polytomy, test all pairwise mergers and choose the one with highest likelihood gain
  let n_resolved = if args.resolve_polytomies {
    resolve_polytomies(graph, partitions)?
  } else {
    0
  };

  if n_resolved > 0 {
    prepare_tree_after_topology_change(graph)?;
    if aln.is_some() {
      run_marginal(graph, partitions, aln)?;
    }
    is_tree_dirty = true;
  }

  // Conditional reconstruction order based on whether tree structure changed:
  // - If changed: update times first (tree -> sequences)
  // - If not changed: update sequences first (sequences -> times)
  let ndiff = if is_tree_dirty {
    run_timetree(graph, partitions).wrap_err_with(|| format!("Timetree inference failed (iteration {i})"))?;

    if aln.is_some() {
      run_marginal(graph, partitions, aln)?;
    }

    // TODO: Implement sequence change tracking
    // What: Count number of ancestral state changes between consecutive iterations
    // Why: Core - needed to detect convergence (ndiff == 0 means no more sequence changes)
    // How: Compare ancestral sequences before/after reconstruction, count differing positions
    0
  } else {
    let sequence_changes = if aln.is_some() {
      run_marginal(graph, partitions, aln)?;
      // TODO: Implement sequence change counting
      // What: Count number of ancestral state changes between consecutive iterations
      // Why: Core - convergence criterion (stop when ndiff == 0)
      // How: Compare ancestral sequences before/after reconstruction, sum differing positions across tree
      0
    } else {
      0
    };

    run_timetree(graph, partitions).wrap_err_with(|| format!("Timetree inference failed (iteration {i})"))?;

    sequence_changes
  };

  Ok((ndiff, n_resolved))
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

pub fn run_timetree_estimation(args: TreetimeTimetreeArgs) -> Result<(), Report> {
  let InputData {
    graph,
    alphabet,
    aln,
    constraints,
  } = load_input_data(&args)?;

  let mut clock_model = infer_clock_model(&args, &graph, &constraints).wrap_err("Failed to infer clock model")?;

  let partitions: Vec<Arc<RwLock<dyn PartitionTreetimeMarginalOps>>> =
    initialize_partitions(&args, &graph, alphabet, aln.as_deref(), &constraints)?;

  run_pre_optimization(&args, &graph, &partitions, aln.as_deref(), &constraints)?;

  run_initial_timetree_inference(&args, &graph, &partitions, aln.as_deref(), &constraints)?;

  let mut optimizer = TimetreeOptimizer::new(args.max_iter, args.tracelog.clone())?;
  while let Some(IterationContext { i }) = optimizer.next_iter() {
    let (ndiff, n_resolved) = run_refinement_iteration(&args, &graph, &partitions, aln.as_deref(), i)?;

    clock_model = update_clock_model(&graph, &constraints, &clock_model)
      .wrap_err_with(|| format!("Failed to update clock model (iteration {i})"))?;

    optimizer
      .record(ndiff, n_resolved, &partitions)
      .wrap_err_with(|| format!("Failed to record convergence metrics (iteration {i})"))?;
  }

  run_post_processing(&args, &graph, &partitions, &constraints, &clock_model)?;

  write_outputs(&args, &graph, &partitions, &constraints, &clock_model)?;

  Ok(())
}

// ============================================================================
// Stub functions for future implementation
// ============================================================================

// ============================================================================
// Rerooting algorithms
// ============================================================================

/// Reroot tree to optimize temporal signal using root-to-tip regression.
///
/// Core: Required for molecular clock analysis when root position affects temporal signal quality.
/// Why: Finding optimal root maximizes correlation between sampling dates and genetic divergence.
/// How: Tests all possible root positions, evaluates regression fit for each.
pub fn reroot_tree(_graph: &GraphAncestral, _constraints: &DateConstraintSet, _method: &str) -> Result<(), Report> {
  todo!("Rerooting: least-squares (minimize residuals), min_dev (minimize variation), oldest (root on oldest sample)")
}

// ============================================================================
// Clock filter and outlier detection
// ============================================================================

/// Detect and mark branches that violate molecular clock assumptions.
///
/// Optional: Improves robustness by excluding outliers from clock model inference.
/// Why: Recombination, selection, or sequencing errors can create non-clock-like branches.
/// How: Root-to-tip regression residuals analyzed via IQD (interquartile distance) threshold.
pub fn clock_filter(_graph: &GraphAncestral, _constraints: &DateConstraintSet, _n_iqd: f64) -> Result<(), Report> {
  todo!("Mark branches as bad_branch=true if |residual| > n_iqd * IQD")
}

/// Identify and report outlier branches that violate molecular clock.
///
/// Core: User needs to know which samples were excluded from analysis.
/// Why: Outliers may indicate data quality issues requiring investigation.
/// How: Compare inferred dates to input constraints, report significant deviations.
pub fn report_bad_branches(_graph: &GraphAncestral, _constraints: &DateConstraintSet) -> Result<(), Report> {
  todo!("Log warnings for branches where |inferred_date - constraint_date| is large")
}

// ============================================================================
// Advanced clock models
// ============================================================================

/// Add coalescent model as population genetic prior on node times.
///
/// Optional: Improves time estimates by incorporating population dynamics.
/// Why: Coalescent theory provides realistic expectations for divergence times.
/// How: Kingman coalescent with exponential waiting times between mergers.
pub fn add_coalescent_model(
  _graph: &GraphAncestral,
  _partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
  _params: &str,
) -> Result<(), Report> {
  todo!("Add merger rate contribution to node time distributions")
}

/// Apply relaxed molecular clock allowing branch-specific rate variation.
///
/// Optional: Accounts for heterogeneity in evolutionary rates across tree.
/// Why: Strict clock assumption may be violated in real data.
/// How: Autocorrelated rates (neighboring branches have similar rates) with penalty terms.
pub fn apply_relaxed_clock(
  _graph: &GraphAncestral,
  _partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
  _params: &[f64],
) -> Result<(), Report> {
  todo!("Optimize branch-specific gamma (rate multiplier) with slack/coupling penalties")
}

// ============================================================================
// Polytomy resolution
// ============================================================================

/// Resolve multifurcations using temporal constraints and sequence data.
///
/// Optional: Improves tree topology when input tree has polytomies.
/// Why: Binary trees are more informative and required for some analyses.
/// How: Greedy optimization of likelihood gain from inserting internal nodes.
pub fn resolve_polytomies(
  _graph: &GraphAncestral,
  _partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
) -> Result<usize, Report> {
  todo!(
    "For each polytomy: test pairwise mergers, choose highest likelihood gain. Returns count of resolved polytomies."
  )
}

/// Reset internal state after tree topology changes.
///
/// Core: Required after polytomy resolution to maintain consistency.
/// Why: Topology changes invalidate cached likelihood values and node attributes.
/// How: Clear cached data, recalculate tree layout, reset branch flags.
pub fn prepare_tree_after_topology_change(_graph: &GraphAncestral) -> Result<(), Report> {
  todo!("Clear cached likelihoods, recalculate dist2root, reset bad_branch flags")
}

// ============================================================================
// Convergence tracking
// ============================================================================

/// Calculate total sequence likelihood across all partitions.
///
/// Core: Component of convergence monitoring.
/// Why: Track changes in sequence reconstruction quality across iterations.
/// How: Sum sequence likelihood from all partitions (joint or marginal depending on mode).
pub fn calculate_sequence_likelihood(_partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>]) -> Option<f64> {
  // TODO: Extract sequence_joint_LH or sequence_marginal_LH from partitions
  // For now return None until partition likelihood tracking is implemented
  None
}

/// Calculate total positional (temporal) likelihood across all partitions.
///
/// Core: Component of convergence monitoring.
/// Why: Track changes in node time inference quality across iterations.
/// How: Sum positional likelihood from all partitions.
pub fn calculate_positional_likelihood(_partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>]) -> Option<f64> {
  // TODO: Extract positional_LH from partitions
  // This represents the likelihood of the inferred node times given temporal constraints
  None
}

/// Calculate total coalescent likelihood across all partitions.
///
/// Optional: Only relevant when coalescent model is active.
/// Why: Track coalescent prior contribution to convergence.
/// How: Sum coalescent likelihood from merger model if present.
pub fn calculate_coalescent_likelihood(_partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>]) -> Option<f64> {
  // TODO: Extract coalescent_joint_LH from merger model if active
  // Returns None when coalescent model not used
  None
}

/// Calculate total likelihood (sum of all components).
///
/// Core: Primary convergence metric.
/// Why: Combined likelihood should increase or stabilize as tree converges.
/// How: Sum sequence, positional, and coalescent likelihoods.
pub fn calculate_total_likelihood(partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>]) -> Option<f64> {
  let seq_lh = calculate_sequence_likelihood(partitions);
  let pos_lh = calculate_positional_likelihood(partitions);
  let coal_lh = calculate_coalescent_likelihood(partitions);

  match (seq_lh, pos_lh, coal_lh) {
    (Some(s), Some(p), Some(c)) => Some(s + p + c),
    (Some(s), Some(p), None) => Some(s + p),
    (Some(s), None, None) => Some(s),
    (None, Some(p), None) => Some(p),
    _ => None,
  }
}

/// Count ancestral sequence changes between iterations.
///
/// Core: Primary convergence criterion (stop when ndiff == 0).
/// Why: No sequence changes indicates reconstruction has stabilized.
/// How: Compare current and previous ancestral sequences, count differing positions.
pub fn count_sequence_changes(
  _previous_states: &AncestralStateSnapshot,
  _current_states: &AncestralStateSnapshot,
) -> usize {
  // TODO: Compare sequences position-by-position across all internal nodes
  // Return total count of positions where ancestral state changed
  0
}

/// Capture current ancestral sequence states for change tracking.
///
/// Core: Required for ndiff calculation between iterations.
/// Why: Need snapshot to compare before/after reconstruction.
/// How: Extract reconstructed sequences from all internal nodes.
pub fn capture_ancestral_states(
  _graph: &GraphAncestral,
  _partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
) -> AncestralStateSnapshot {
  // TODO: Extract ancestral sequences from partition node data
  // Map node keys to their current sequence states
  BTreeMap::new()
}

// ============================================================================
// Sensitivity analysis
// ============================================================================

/// Assess sensitivity of timetree to clock rate uncertainty.
///
/// Optional: Quantifies robustness of time estimates.
/// Why: Clock rate has confidence intervals; propagating uncertainty improves reliability.
/// How: Re-run timetree with rate ± std_dev, compare resulting node times.
pub fn calc_rate_susceptibility(
  _graph: &GraphAncestral,
  _partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
  _constraints: &DateConstraintSet,
  _clock_model: &ClockModel,
) -> Result<(), Report> {
  todo!("Run timetree with rate±σ, store alternative time estimates")
}

/// Extract confidence intervals from marginal posterior distributions.
///
/// Core: Quantifies uncertainty in inferred node times.
/// Why: Point estimates alone don't convey reliability of time inference.
/// How: Compute HPD (highest posterior density) or quantiles from marginal distributions.
pub fn extract_confidence_intervals(
  _graph: &GraphAncestral,
  _partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
) -> Result<(), Report> {
  todo!("For each node: compute 95% HPD interval from marginal_pos_LH distribution")
}

// ============================================================================
// Output functions
// ============================================================================

/// Write inferred node dates to file.
///
/// Core: Primary output of timetree inference.
/// Why: Users need calendar dates for each node.
/// How: TSV file with node_name, numdate, time_before_present.
pub fn write_node_dates(
  _graph: &GraphAncestral,
  _partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
  _out_base: &Path,
) -> Result<(), Report> {
  todo!("Write node dates to TSV file")
}

/// Write confidence intervals for node dates.
///
/// Core: Essential for assessing reliability of time estimates.
/// Why: Quantifies uncertainty in molecular dating.
/// How: TSV file with node_name, lower_bound, median, upper_bound.
pub fn write_confidence_intervals(
  _graph: &GraphAncestral,
  _partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
  _out_base: &Path,
) -> Result<(), Report> {
  todo!("Write confidence intervals to TSV file")
}

/// Write clock model parameters to file.
///
/// Core: Documents inferred molecular clock rate.
/// Why: Clock rate is key parameter for interpreting divergence times.
/// How: JSON file with rate, intercept, R², confidence intervals.
pub fn write_clock_model(_clock_model: &ClockModel, _out_base: &Path) -> Result<(), Report> {
  todo!("Write clock model to JSON file")
}

/// Plot root-to-tip distance vs sampling date regression.
///
/// Optional: Visual quality control for temporal signal.
/// Why: Diagnostic for clock-like evolution and outlier detection.
/// How: Scatter plot with regression line, residuals, R² annotation.
pub fn plot_root_to_tip(
  _graph: &GraphAncestral,
  _constraints: &DateConstraintSet,
  _out_base: &Path,
) -> Result<(), Report> {
  todo!("Generate root-to-tip regression plot")
}

/// Plot time-scaled phylogenetic tree.
///
/// Optional: Visual representation of inferred divergence times.
/// Why: Tree topology combined with temporal axis aids interpretation.
/// How: Horizontal layout with x-axis as calendar time, branches scaled by divergence times.
pub fn plot_time_tree(_graph: &GraphAncestral, _out_base: &Path) -> Result<(), Report> {
  todo!("Generate time-scaled tree visualization")
}
