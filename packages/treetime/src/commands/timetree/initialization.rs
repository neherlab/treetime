use crate::representation::edge_timetree::EdgeTimetree;
use crate::representation::node_timetree::NodeTimetree;
use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::get_common_length;
use crate::representation::infer_dense::infer_dense;
use crate::commands::timetree::args::{BranchLengthMode, TreetimeTimetreeArgs};
use crate::commands::timetree::data::date_constraints::DateConstraintSet;
use crate::commands::timetree::data::date_constraints::load_date_constraints;
use crate::gtr::get_gtr::{JC69Params, jc69};
use crate::io::fasta::{FastaRecord, read_many_fasta};
use crate::io::nwk::nwk_read_file;
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::partition_timetree::{GraphTimetree, PartitionTreetimeMarginalOps, PartitionTreetimeMarginalVec};
use crate::representation::partition_timetree_dense::PartitionTimetreeDense;
use crate::representation::partition_timetree_sparse::PartitionTimetreeSparse;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use maplit::btreemap;
use parking_lot::RwLock;
use std::sync::Arc;

pub struct InputData {
  pub graph: GraphTimetree,
  pub alphabet: Alphabet,
  pub aln: Option<Vec<FastaRecord>>,
  pub constraints: DateConstraintSet,
}

pub fn load_input_data(args: &TreetimeTimetreeArgs) -> Result<InputData, Report> {
  let graph_ancestral: GraphAncestral = if let Some(tree_path) = &args.tree {
    nwk_read_file(tree_path).wrap_err("Failed to load tree from file")?
  } else {
    todo!("Tree inference from alignment not yet implemented")
  };

  // Convert GraphAncestral to GraphTimetree
  let graph: GraphTimetree = graph_ancestral.to_graph_timetree()?;

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

pub fn initialize_partitions(
  args: &TreetimeTimetreeArgs,
  _graph: &GraphTimetree,
  alphabet: Alphabet,
  aln: Option<&[FastaRecord]>,
  _constraints: &DateConstraintSet,
) -> Result<PartitionTreetimeMarginalVec, Report> {
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
    })) as Arc<RwLock<dyn PartitionTreetimeMarginalOps<NodeTimetree, EdgeTimetree>>>]
  } else {
    [Arc::new(RwLock::new(PartitionTimetreeDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      sequence_length,
      nodes: btreemap! {},
      edges: btreemap! {},
    })) as Arc<RwLock<dyn PartitionTreetimeMarginalOps<NodeTimetree, EdgeTimetree>>>]
  }
  .into_iter()
  .collect_vec();

  Ok(partitions)
}
