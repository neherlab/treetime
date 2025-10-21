use crate::commands::ancestral::marginal_unified::run_marginal;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::inference::runner::run_timetree;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::io::fasta::FastaRecord;
use crate::representation::graph_ancestral::GraphAncestral;
use eyre::{Report, WrapErr};
use parking_lot::RwLock;
use std::sync::Arc;

#[allow(clippy::useless_let_if_seq)]
pub fn run_refinement_iteration(
  args: &TreetimeTimetreeArgs,
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll>>],
  aln: Option<&[FastaRecord]>,
  i: usize,
) -> Result<(usize, usize), Report> {
  let mut is_tree_dirty = false;

  if let Some(_coalescent_params) = &args.coalescent {
    todo!("add_coalescent_model not yet implemented");
  }

  if !args.relax.is_empty() {
    todo!("apply_relaxed_clock not yet implemented");
  }

  let n_resolved = if args.resolve_polytomies {
    todo!("resolve_polytomies not yet implemented");
  } else {
    0
  };

  if n_resolved > 0 {
    is_tree_dirty = true;
  }

  if is_tree_dirty {
    run_timetree(graph, partitions).wrap_err_with(|| format!("Timetree inference failed (iteration {i})"))?;

    if aln.is_some() {
      run_marginal(graph, partitions, aln)?;
    }
  } else {
    if aln.is_some() {
      run_marginal(graph, partitions, aln)?;
    }

    run_timetree(graph, partitions).wrap_err_with(|| format!("Timetree inference failed (iteration {i})"))?;
  }

  let ndiff = 0;

  Ok((ndiff, n_resolved))
}
