use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::get_common_length;
use crate::commands::timetree::args::{BranchLengthMode, TreetimeTimetreeArgs};
use crate::commands::timetree::data::date_constraints::load_date_constraints;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::gtr::get_gtr::{JC69Params, jc69};
use crate::io::dates_csv::read_dates;
use crate::io::fasta::{FastaRecord, read_many_fasta};
use crate::io::nwk::nwk_read_file;
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::infer_dense::infer_dense;
use crate::representation::partition_marginal_dense::PartitionMarginalDense;
use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
use eyre::{Report, WrapErr};
use maplit::btreemap;
use parking_lot::RwLock;
use std::sync::Arc;

pub struct InputData {
  pub graph: GraphAncestral,
  pub alphabet: Alphabet,
  pub aln: Option<Vec<FastaRecord>>,
}

pub fn load_input_data(args: &TreetimeTimetreeArgs) -> Result<InputData, Report> {
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

  if let Some(dates_path) = &args.dates {
    let dates = read_dates(dates_path, &args.name_column, &args.date_column).wrap_err("When reading dates")?;
    load_date_constraints(&dates, &graph).wrap_err("Failed to load date constraints")?;
  }

  Ok(InputData { graph, alphabet, aln })
}

pub fn initialize_partitions(
  args: &TreetimeTimetreeArgs,
  _graph: &GraphAncestral,
  alphabet: Alphabet,
  aln: Option<&[FastaRecord]>,
) -> Result<Vec<Arc<RwLock<dyn PartitionTimetreeAll>>>, Report> {
  let dense = args.dense.unwrap_or_else(infer_dense);
  let length = if let Some(aln_data) = aln {
    get_common_length(aln_data)?
  } else {
    args
      .sequence_length
      .ok_or_else(|| eyre::eyre!("sequence_length required when no alignment provided"))?
  };

  if !dense {
    let partition: Arc<RwLock<dyn PartitionTimetreeAll>> =
      Arc::new(RwLock::new(PartitionMarginalSparse {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet,
        length,
        nodes: btreemap! {},
        edges: btreemap! {},
      }));
    Ok(vec![partition])
  } else {
    let partition: Arc<RwLock<dyn PartitionTimetreeAll>> =
      Arc::new(RwLock::new(PartitionMarginalDense {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet,
        length,
        nodes: btreemap! {},
        edges: btreemap! {},
      }));
    Ok(vec![partition])
  }
}
