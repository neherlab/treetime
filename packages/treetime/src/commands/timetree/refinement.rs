use crate::representation::edge_timetree::EdgeTimetree;
use crate::representation::node_timetree::NodeTimetree;
use crate::commands::ancestral::marginal_unified::run_marginal;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::inference::runner::run_timetree;
use crate::commands::timetree::optimization::coalescent::add_coalescent_model;
use crate::commands::timetree::optimization::polytomy::{prepare_tree_after_topology_change, resolve_polytomies};
use crate::commands::timetree::optimization::relaxed_clock::apply_relaxed_clock;
use crate::io::fasta::FastaRecord;
use crate::representation::partition_timetree::{GraphTimetree, PartitionTreetimeMarginalOps};
use eyre::{Report, WrapErr};
use parking_lot::RwLock;
use std::sync::Arc;

#[allow(clippy::useless_let_if_seq)]
pub fn run_refinement_iteration(
  args: &TreetimeTimetreeArgs,
  graph: &GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps<NodeTimetree, EdgeTimetree>>>],
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
  if is_tree_dirty {
    run_timetree(graph, partitions).wrap_err_with(|| format!("Timetree inference failed (iteration {i})"))?;

    if aln.is_some() {
      run_marginal(graph, partitions, aln)?;
    }

    // TODO: Implement sequence change tracking
    // What: Count number of ancestral state changes between consecutive iterations
    // Why: Core - needed to detect convergence (ndiff == 0 means no more sequence changes)
    // How: Compare ancestral sequences before/after reconstruction, count differing positions
  } else {
    if aln.is_some() {
      run_marginal(graph, partitions, aln)?;
      // TODO: Implement sequence change counting
      // What: Count number of ancestral state changes between consecutive iterations
      // Why: Core - convergence criterion (stop when ndiff == 0)
      // How: Compare ancestral sequences before/after reconstruction, sum differing positions across tree
    }

    run_timetree(graph, partitions).wrap_err_with(|| format!("Timetree inference failed (iteration {i})"))?;
  }

  let ndiff = 0;

  Ok((ndiff, n_resolved))
}
