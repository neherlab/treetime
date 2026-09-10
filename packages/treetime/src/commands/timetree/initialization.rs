use crate::alphabet::alphabet::Alphabet;
use crate::clock::date_constraints::load_date_constraints;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::make_error;
use crate::make_report;
use crate::optimize::params::BranchLengthMode;
use crate::partition::timetree::partition::GraphTimetree;
use crate::seq::gap_fill::apply_gap_fill;
use eyre::{Report, WrapErr};
use treetime_graph::value_maps::node_names;
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
  let names = node_names(&graph);
  let input_leaf_order = graph
    .get_leaves()
    .into_iter()
    .map(|leaf| {
      let key = leaf.read_arc().key();
      names[&key]
        .clone()
        .ok_or_else(|| make_report!("Leaf node {key} has no name"))
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
      &args.metadata_id.metadata_delimiters,
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

  Ok(InputData {
    graph,
    input_leaf_order,
    alphabet,
    aln,
    dates,
  })
}
