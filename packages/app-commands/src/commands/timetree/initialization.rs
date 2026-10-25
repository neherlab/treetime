use crate::commands::shared::dates_input::read_input_dates;
use crate::commands::shared::leaf_order::leaf_order;
use crate::commands::shared::sequence_inputs::{AlignmentSource, SequenceSources};
use crate::commands::shared::tree_input::read_input_tree;
use crate::commands::timetree::args::{TreetimeTimetreeArgs, check_branch_length_settings};
use crate::commands::timetree::branch_length_mode::{ResolvedBranchLengthMode, resolve_branch_length_mode};
use eyre::{Report, WrapErr};
use treetime::partition::resolve::PartitionPlan;
use treetime::progress::LogSink;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_io::dates_csv::DateConstraint;

pub(crate) fn load_input_data(args: &TreetimeTimetreeArgs, log: &dyn LogSink) -> Result<InputData, Report> {
  let nwk_parsed = read_input_tree(&args.tree, args.newick_dialect, log).wrap_err("Failed to load tree from file")?;
  let names = nwk_parsed.names;
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let input_leaf_order = leaf_order(&graph);
  let has_alignment = !args.sequences.alignments.is_empty();
  let branch_length_mode = resolve_branch_length_mode(args.branch_lengths.mode, has_alignment, &branch_lengths, log);
  check_branch_length_settings(&args.marginal_only_settings, &branch_length_mode)?;

  let SequenceSources {
    alignments: sequences,
    partitions,
  } = args
    .sequences
    .read_partitions(&args.sequence_model.model, &graph, &names, log)?;

  let dates = read_input_dates(&args.dates, &graph, &names, log)?;

  Ok(InputData {
    graph,
    names,
    branch_lengths,
    branch_length_mode,
    input_leaf_order,
    sequences,
    partitions,
    dates,
  })
}

pub(crate) struct InputData {
  pub graph: Graph,
  pub names: NodeNames,
  pub branch_lengths: EdgeMap<Option<f64>>,
  pub branch_length_mode: ResolvedBranchLengthMode,
  pub input_leaf_order: Vec<GraphNodeKey>,
  pub sequences: Vec<AlignmentSource>,
  pub partitions: Vec<PartitionPlan>,
  pub dates: NodeMap<DateConstraint>,
}
