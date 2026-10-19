use crate::commands::shared::dates_input::read_input_dates;
use crate::commands::shared::leaf_order::leaf_order;
use crate::commands::shared::sequence_inputs::{AlignmentSource, SourceSettings, read_sequence_sources};
use crate::commands::shared::tree_input::read_input_tree;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::attach::MissingLeafLimit;
use treetime::progress::LogSink;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::dates_csv::DateConstraint;

pub(crate) fn load_input_data(args: &TreetimeTimetreeArgs, log: &dyn LogSink) -> Result<InputData, Report> {
  let nwk_parsed =
    read_input_tree(&args.tree, args.tree_dialect.dialect(), log).wrap_err("Failed to load tree from file")?;
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let input_leaf_order = leaf_order(&graph);

  let alphabet = Alphabet::new(args.alphabet_args.alphabet_name().unwrap_or_default())?;

  let sequences = if args.alignment.alignment.is_empty() {
    vec![]
  } else {
    let settings = SourceSettings {
      alphabet: &alphabet,
      gap_fill: Some(args.gap_fill_args.effective_gap_fill()),
      missing_leaves: MissingLeafLimit::Enforced,
    };
    read_sequence_sources(&args.alignment, &settings, &graph, &names, log)?
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
    sequences,
    dates,
  })
}

pub(crate) struct InputData {
  pub graph: Graph,
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub input_leaf_order: Vec<GraphNodeKey>,
  pub sequences: Vec<AlignmentSource>,
  pub dates: Option<BTreeMap<GraphNodeKey, DateConstraint>>,
}
