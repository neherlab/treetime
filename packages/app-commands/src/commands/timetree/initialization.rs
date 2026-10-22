use crate::commands::shared::dates_input::read_input_dates;
use crate::commands::shared::leaf_order::leaf_order;
use crate::commands::shared::sequence_inputs::{AlignmentSource, SourceSettings, read_sequence_sources};
use crate::commands::shared::tree_input::read_input_tree;
use crate::commands::timetree::args::{TreetimeTimetreeArgs, check_branch_length_settings};
use crate::commands::timetree::branch_length_mode::{ResolvedBranchLengthMode, resolve_branch_length_mode};
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use treetime::alphabet::alphabet::Alphabet;
use treetime::progress::LogSink;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_io::dates_csv::DateConstraint;

pub(crate) fn load_input_data(args: &TreetimeTimetreeArgs, log: &dyn LogSink) -> Result<InputData, Report> {
  let nwk_parsed =
    read_input_tree(&args.tree, args.tree_dialect.dialect(), log).wrap_err("Failed to load tree from file")?;
  let names = nwk_parsed.names;
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let input_leaf_order = leaf_order(&graph);
  let has_alignment = !args.alignment.alignment.is_empty();
  let branch_length_mode = resolve_branch_length_mode(args.branch_length_mode, has_alignment, &branch_lengths, log);
  check_branch_length_settings(args, &branch_length_mode)?;

  let alphabet = Alphabet::new(args.alphabet_args.alphabet_name().unwrap_or_default())?;

  let sequences = if has_alignment {
    let settings = SourceSettings {
      alphabet: &alphabet,
      gap_fill: Some(args.gap_fill_args.gap_fill()),
      missing_leaves: args.missing_sequences.limit(),
    };
    read_sequence_sources(&args.alignment, &settings, &graph, &names, log)?
  } else {
    vec![]
  };

  let dates = args
    .metadata
    .as_deref()
    .map(|path| {
      read_input_dates(
        path,
        &args.metadata_id,
        args.date_column.date_column.as_deref(),
        &args.date_column.date_format,
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
    branch_length_mode,
    input_leaf_order,
    sequences,
    dates,
  })
}

pub(crate) struct InputData {
  pub graph: Graph,
  pub names: NodeNames,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub branch_length_mode: ResolvedBranchLengthMode,
  pub input_leaf_order: Vec<GraphNodeKey>,
  pub sequences: Vec<AlignmentSource>,
  pub dates: Option<BTreeMap<GraphNodeKey, DateConstraint>>,
}
