use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::get_common_length;
use crate::commands::ancestral::marginal::update_marginal;
use crate::commands::clock::date_constraints::load_date_constraints;
use crate::commands::timetree::args::{BranchLengthMode, TreetimeTimetreeArgs};
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::commands::timetree::utils::initialize_node_divergences;
use crate::gtr::get_gtr::{
  GtrModelName, JC69Params, get_gtr_by_name, get_gtr_dense, get_gtr_sparse, jc69, log_gtr, write_gtr_json,
};
use crate::make_error;
use crate::make_report;
use crate::representation::algo::infer_dense::infer_dense;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition::timetree::{GraphTimetree, PartitionTimetreeAllVec};
use crate::representation::partition::traits::PartitionMarginalOps;
use crate::representation::payload::timetree::EdgeTimetree;
use crate::representation::payload::timetree::NodeTimetree;
use eyre::{Report, WrapErr};
use log::info;
use maplit::btreemap;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_io::dates_csv::read_dates;
use treetime_io::fasta::{FastaRecord, read_many_fasta};
use treetime_io::nwk::nwk_read_file;
use treetime_primitives::seq;

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
  let model_name = args.gtr;
  let length = if let Some(aln_data) = aln {
    get_common_length(aln_data)?
  } else {
    args
      .sequence_length
      .ok_or_else(|| make_report!("sequence_length required when no alignment provided"))?
  };

  // For named models, construct the specified GTR directly.
  // For Infer, start with JC69 placeholder, then infer the real model from
  // reconstructed sequences (sparse: after Fitch compression, dense: after
  // an initial marginal pass). Matches v0 behavior where GTR inference
  // happens before rerooting.
  let initial_gtr = if model_name == GtrModelName::Infer {
    info!("GTR model: infer (starting with JC69 placeholder)");
    let gtr = jc69(JC69Params::default())?;
    log_gtr(&gtr, GtrModelName::JC69);
    gtr
  } else {
    info!("GTR model: {model_name}");
    let gtr = get_gtr_by_name(model_name).wrap_err_with(|| format!("When creating GTR model '{model_name}'"))?;
    log_gtr(&gtr, model_name);
    gtr
  };

  if !dense {
    let aln_data = aln.ok_or_else(|| make_report!("Alignment required for sparse marginal reconstruction"))?;

    let sparse_partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: initial_gtr,
      alphabet,
      length,
      root_sequence: seq![],
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    crate::commands::ancestral::fitch::compress_sequences(graph, std::slice::from_ref(&sparse_partition), aln_data)?;
    sparse_partition.write_arc().extract_root_sequence(graph);

    // For Infer: Fitch compression populated mutation counts, infer real GTR
    if model_name == GtrModelName::Infer {
      let gtr = get_gtr_sparse(&model_name, &sparse_partition, graph)
        .wrap_err("When inferring GTR model from sparse partition")?;
      sparse_partition.write_arc().gtr = gtr;
    }

    write_gtr_json(&sparse_partition.read_arc().gtr, model_name, &args.outdir, None)?;

    let partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> = sparse_partition;
    Ok(vec![partition])
  } else {
    let dense_partition = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: initial_gtr,
      alphabet,
      length,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    // For Infer: run an initial marginal pass to populate profiles, then
    // infer the real GTR from the branch joint distributions. The profiles
    // will be recomputed later by initialize_marginal with the inferred GTR.
    if model_name == GtrModelName::Infer {
      let aln_data = aln.ok_or_else(|| make_report!("Alignment required for dense GTR inference"))?;
      dense_partition.write_arc().attach_sequences(graph, aln_data)?;
      update_marginal(graph, std::slice::from_ref(&dense_partition))?;
      let gtr = get_gtr_dense(&model_name, &dense_partition, graph)
        .wrap_err("When inferring GTR model from dense partition")?;
      dense_partition.write_arc().gtr = gtr;
    }

    write_gtr_json(&dense_partition.read_arc().gtr, model_name, &args.outdir, None)?;

    let partition: Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>> = dense_partition;
    Ok(vec![partition])
  }
}
