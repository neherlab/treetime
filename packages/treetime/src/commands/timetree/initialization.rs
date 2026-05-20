use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch::create_fitch_partition;
use crate::ancestral::gtr_inference::infer_gtr_fitch;
use crate::clock::date_constraints::load_date_constraints;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::gtr::get_gtr::{GtrModelName, get_gtr_by_name, log_gtr, write_gtr_json};
use crate::make_error;
use crate::make_report;
use crate::optimize::params::BranchLengthMode;
use crate::partition::algo::infer_dense::infer_dense;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::timetree::{GraphTimetree, PartitionTimetreeAllVec};
use crate::partition::traits::PartitionTimetreeAll;
use crate::payload::timetree::EdgeTimetree;
use crate::payload::timetree::NodeTimetree;
use crate::seq::alignment::get_common_length;
use crate::seq::gap_fill::apply_gap_fill;
use crate::timetree::utils::initialize_node_divergences;
use eyre::{Report, WrapErr};
use log::info;
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

  let alphabet = Alphabet::new(args.alphabet)?;

  // Load alignment sequences (optional if using input branch lengths only)
  let aln = if !args.input_fastas.is_empty() {
    let mut records = read_many_fasta(&args.input_fastas, &alphabet)?;
    let gap_fill_mode = args.effective_gap_fill();
    for record in &mut records {
      apply_gap_fill(&mut record.seq, gap_fill_mode, alphabet.gap(), alphabet.unknown());
    }
    Some(records)
  } else if args.branch_length_mode != BranchLengthMode::Input {
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
  let model_name = args.gtr;
  let length = if let Some(aln_data) = aln {
    get_common_length(aln_data)?
  } else {
    args
      .sequence_length
      .ok_or_else(|| make_report!("sequence_length required when no alignment provided"))?
  };

  if !dense {
    let aln_data = aln.ok_or_else(|| make_report!("Alignment required for sparse marginal reconstruction"))?;

    let fitch = create_fitch_partition(graph, 0, alphabet, aln_data)?;
    let gtr = match model_name {
      GtrModelName::Infer => infer_gtr_fitch(&fitch, graph)?,
      _ => get_gtr_by_name(model_name)?,
    };
    log_gtr(&gtr, model_name);
    let partition = fitch.into_marginal_sparse(gtr, graph)?;

    write_gtr_json(&partition.gtr, model_name, &args.outdir, None)?;

    let sparse_partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> =
      Arc::new(RwLock::new(partition));
    Ok(vec![sparse_partition])
  } else if model_name == GtrModelName::Infer {
    let aln_data = aln.ok_or_else(|| make_report!("Alignment required for dense GTR inference"))?;
    let fitch = create_fitch_partition(graph, 0, alphabet, aln_data)?;
    let gtr = infer_gtr_fitch(&fitch, graph)?;
    log_gtr(&gtr, model_name);
    let partition = fitch.into_marginal_dense(gtr);

    write_gtr_json(&partition.data.gtr, model_name, &args.outdir, None)?;

    let dense_partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> =
      Arc::new(RwLock::new(partition));
    Ok(vec![dense_partition])
  } else {
    info!("GTR model: {model_name}");
    let gtr = get_gtr_by_name(model_name).wrap_err_with(|| format!("When creating GTR model '{model_name}'"))?;
    log_gtr(&gtr, model_name);
    let partition = PartitionMarginalDense::new(0, gtr, alphabet, length);

    write_gtr_json(&partition.data.gtr, model_name, &args.outdir, None)?;

    let dense_partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> =
      Arc::new(RwLock::new(partition));
    Ok(vec![dense_partition])
  }
}
