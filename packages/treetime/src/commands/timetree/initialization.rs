use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::get_common_length;
use crate::commands::clock::date_constraints::load_date_constraints;
use crate::commands::timetree::args::{BranchLengthMode, TreetimeTimetreeArgs};
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::commands::timetree::utils::initialize_node_divergences;
use crate::gtr::get_gtr::{GtrModelName, JC69Params, jc69, log_gtr};
use crate::make_error;
use crate::make_report;
use crate::representation::algo::infer_dense::infer_dense;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition::timetree::{GraphTimetree, PartitionTimetreeAllVec};
use crate::representation::payload::timetree::EdgeTimetree;
use crate::representation::payload::timetree::NodeTimetree;
use eyre::{Report, WrapErr};
use maplit::btreemap;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_io::dates_csv::read_dates;
use treetime_io::fasta::{FastaRecord, read_many_fasta};
use treetime_io::nwk::nwk_read_file;

pub struct InputData {
  pub graph: GraphTimetree,
  pub alphabet: Alphabet,
  pub aln: Option<Vec<FastaRecord>>,
}

pub fn load_input_data(args: &TreetimeTimetreeArgs) -> Result<InputData, Report> {
  let graph: GraphTimetree = if let Some(tree_path) = &args.tree {
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
    return make_error!(
      "Alignment required when branch_length_mode is not 'input'. \
       Provide FASTA files or use --branch-length-mode=input"
    );
  } else {
    None
  };

  if let Some(dates_path) = &args.dates {
    let dates = read_dates(dates_path, &args.name_column, &args.date_column).wrap_err("When reading dates")?;
    load_date_constraints(&dates, &graph).wrap_err("Failed to load date constraints")?;
  }

  // Calculate divergence distances from root to all nodes
  initialize_node_divergences(&graph);

  Ok(InputData { graph, alphabet, aln })
}

pub fn initialize_partitions(
  args: &TreetimeTimetreeArgs,
  graph: &GraphTimetree,
  alphabet: Alphabet,
  aln: Option<&[FastaRecord]>,
) -> Result<PartitionTimetreeAllVec, Report> {
  let dense = args.dense.unwrap_or_else(infer_dense);
  let length = if let Some(aln_data) = aln {
    get_common_length(aln_data)?
  } else {
    args
      .sequence_length
      .ok_or_else(|| make_report!("sequence_length required when no alignment provided"))?
  };

  if !dense {
    let aln_data = aln.ok_or_else(|| make_report!("Alignment required for sparse marginal reconstruction"))?;

    let gtr = jc69(JC69Params::default())?;
    log_gtr(&gtr, GtrModelName::JC69);

    let sparse_partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr,
      alphabet,
      length,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    crate::commands::ancestral::fitch::compress_sequences(graph, std::slice::from_ref(&sparse_partition), aln_data)?;

    let partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> = sparse_partition;
    Ok(vec![partition])
  } else {
    let gtr = jc69(JC69Params::default())?;
    log_gtr(&gtr, GtrModelName::JC69);

    let partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> =
      Arc::new(RwLock::new(PartitionMarginalDense {
        index: 0,
        gtr,
        alphabet,
        length,
        nodes: btreemap! {},
        edges: btreemap! {},
      }));
    Ok(vec![partition])
  }
}
