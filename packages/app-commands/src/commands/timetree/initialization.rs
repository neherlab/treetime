use crate::commands::shared::alignment::read_alignment;
use crate::commands::shared::leaf_order::leaf_order;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use treetime::alphabet::alphabet::Alphabet;
use treetime::clock::date_constraints::load_date_constraints;
use treetime::make_error;
use treetime::optimize::params::BranchLengthMode;
use treetime::progress::LogSink;
use treetime::seq::gap_fill::apply_gap_fill;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::dates_csv::{DatesMap, metadata_read_file};
use treetime_io::fasta::FastaRecord;
use treetime_io::nwk::nwk_read_file;

pub(crate) fn load_input_data(args: &TreetimeTimetreeArgs, log: &dyn LogSink) -> Result<InputData, Report> {
  let nwk_parsed = nwk_read_file(&args.tree).wrap_err("Failed to load tree from file")?;
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let input_leaf_order = leaf_order(&graph, &names)?;

  let alphabet = Alphabet::new(args.alphabet_args.alphabet_name().unwrap_or_default())?;

  let aln = if !args.alignment.alignment.is_empty() {
    let mut records = read_alignment(&args.alignment.alignment, &alphabet)?;
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
    let dates = metadata_read_file(
      dates_path,
      &args.metadata_id.metadata_delimiters,
      &args.metadata_id.metadata_id_columns,
      None,
      args.date_column.date_column.as_deref(),
    )
    .and_then(|table| table.dates())
    .wrap_err("When reading dates")?;
    load_date_constraints(&dates, &graph, &names, log).wrap_err("Failed to load date constraints")?;
    Some(dates)
  } else {
    None
  };

  Ok(InputData {
    graph,
    names,
    branch_lengths,
    input_leaf_order,
    alphabet,
    aln,
    dates,
  })
}

pub(crate) struct InputData {
  pub graph: Graph,
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub input_leaf_order: Vec<String>,
  pub alphabet: Alphabet,
  pub aln: Option<Vec<FastaRecord>>,
  pub dates: Option<DatesMap>,
}
