use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch::create_fitch_partition;
use crate::ancestral::gtr_inference::infer_gtr_fitch;
use crate::clock::date_constraints::load_date_constraints;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::gtr::get_gtr::{GtrModelName, get_gtr_by_name, log_gtr};
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
use treetime_graph::node::Named;
use treetime_io::dates_csv::{DatesMap, read_dates};
use treetime_io::fasta::{FastaRecord, read_many_fasta};
use treetime_io::nwk::nwk_read_file;

pub struct InputData {
  pub graph: GraphTimetree,
  pub input_leaf_order: Vec<String>,
  pub alphabet: Alphabet,
  pub aln: Option<Vec<FastaRecord>>,
  /// Parsed date constraints, retained for node data JSON output (`raw_date`,
  /// `date_inferred`). `None` when no dates file was provided.
  pub dates: Option<DatesMap>,
}

pub fn load_input_data(args: &TreetimeTimetreeArgs) -> Result<InputData, Report> {
  let graph: GraphTimetree = if let Some(tree_path) = &args.tree {
    nwk_read_file(tree_path).wrap_err("Failed to load tree from file")?
  } else {
    todo!("Tree inference from alignment not yet implemented")
  };
  let input_leaf_order = graph
    .get_leaves()
    .into_iter()
    .map(|leaf| {
      let leaf = leaf.read_arc();
      leaf
        .payload()
        .read_arc()
        .name()
        .map(|name| name.as_ref().to_owned())
        .ok_or_else(|| make_report!("Leaf node {} has no name", leaf.key()))
    })
    .collect::<Result<Vec<_>, _>>()?;

  let alphabet = Alphabet::new(args.alphabet_args.alphabet.unwrap_or_default())?;

  let aln = if !args.alignment.alignment.is_empty() {
    let mut records = read_many_fasta(&args.alignment.alignment, &alphabet)?;
    let gap_fill_mode = args.gap_fill_args.effective_gap_fill();
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

  let dates = if let Some(dates_path) = &args.metadata {
    let dates = read_dates(
      dates_path,
      &args.metadata_id.metadata_id_columns,
      &None,
      &args.date_column_args.date_column,
    )
    .wrap_err("When reading dates")?;
    load_date_constraints(&dates, &graph).wrap_err("Failed to load date constraints")?;
    Some(dates)
  } else {
    None
  };

  // Calculate divergence distances from root to all nodes
  initialize_node_divergences(&graph)?;

  Ok(InputData {
    graph,
    input_leaf_order,
    alphabet,
    aln,
    dates,
  })
}

pub fn initialize_partitions(
  args: &TreetimeTimetreeArgs,
  graph: &GraphTimetree,
  alphabet: Alphabet,
  aln: Option<&[FastaRecord]>,
) -> Result<PartitionTimetreeAllVec, Report> {
  let dense = args.dense.unwrap_or_else(infer_dense);
  let model_name = args.model_args.model;
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

    let sparse_partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> =
      Arc::new(RwLock::new(partition));
    Ok(vec![sparse_partition])
  } else if model_name == GtrModelName::Infer {
    let aln_data = aln.ok_or_else(|| make_report!("Alignment required for dense GTR inference"))?;
    let fitch = create_fitch_partition(graph, 0, alphabet, aln_data)?;
    let gtr = infer_gtr_fitch(&fitch, graph)?;
    log_gtr(&gtr, model_name);
    let partition = fitch.into_marginal_dense(gtr);

    let dense_partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> =
      Arc::new(RwLock::new(partition));
    Ok(vec![dense_partition])
  } else {
    info!("GTR model: {model_name}");
    let gtr = get_gtr_by_name(model_name).wrap_err_with(|| format!("When creating GTR model '{model_name}'"))?;
    log_gtr(&gtr, model_name);
    let partition = PartitionMarginalDense::new(0, gtr, alphabet, length);

    let dense_partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> =
      Arc::new(RwLock::new(partition));
    Ok(vec![dense_partition])
  }
}
