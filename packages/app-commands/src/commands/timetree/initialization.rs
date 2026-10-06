use crate::commands::shared::alignment::{PairedAlignment, pair_alignment, read_alignment};
use crate::commands::shared::dates_input::read_input_dates;
use crate::commands::shared::leaf_order::leaf_order;
use crate::commands::shared::tree_input::read_input_tree;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use treetime::alphabet::alphabet::Alphabet;
use treetime::make_error;
use treetime::optimize::params::BranchLengthMode;
use treetime::progress::LogSink;
use treetime::seq::gap_fill::apply_gap_fill;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::dates_csv::DateConstraint;

pub(crate) fn load_input_data(args: &TreetimeTimetreeArgs, log: &dyn LogSink) -> Result<InputData, Report> {
  let nwk_parsed = read_input_tree(&args.tree, log).wrap_err("Failed to load tree from file")?;
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let input_leaf_order = leaf_order(&graph);

  let alphabet = Alphabet::new(args.alphabet_args.alphabet_name().unwrap_or_default())?;

  let sequences = if !args.alignment.alignment.is_empty() {
    let mut records = read_alignment(&args.alignment.alignment, &alphabet)?;
    let gap_fill_mode = args.gap_fill_args.effective_gap_fill();
    for record in &mut records {
      apply_gap_fill(&mut record.seq, gap_fill_mode, alphabet.gap(), alphabet.unknown());
    }
    Some(pair_alignment(records, &args.alignment.alignment, &graph, &names, log))
  } else if args.branch_length_mode != BranchLengthMode::Input {
    return make_error!(
      "Alignment required when branch_length_mode is not 'input'. \
       Provide FASTA files or use --branch-length-mode=input"
    );
  } else {
    None
  };

  let dates = args
    .metadata
    .as_deref()
    .map(|path| {
      read_input_dates(
        path,
        &args.metadata_id,
        args.date_column.date_column.as_deref(),
        &graph,
        &names,
        log,
      )
    })
    .transpose()?;

  Ok(InputData {
    graph,
    names,
    branch_lengths,
    input_leaf_order,
    alphabet,
    sequences,
    dates,
  })
}

pub(crate) struct InputData {
  pub graph: Graph,
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub input_leaf_order: Vec<GraphNodeKey>,
  pub alphabet: Alphabet,
  pub sequences: Option<PairedAlignment>,
  pub dates: Option<BTreeMap<GraphNodeKey, DateConstraint>>,
}
