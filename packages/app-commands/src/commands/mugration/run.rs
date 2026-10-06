use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::commands::shared::input_warnings::warn_duplicate_names;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::shared::tree_input::read_input_tree;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeTraits};
use app_output::augur_node_data_traits::write_augur_node_data_traits;
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
use app_output::trait_tables::{write_trait_confidence_csv, write_traits_csv};
use app_output::tree_output::{tree_view_for_outputs, write_graph_outputs, write_tree_outputs};
use eyre::Report;
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use treetime::cancel::Cancel;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput};
use treetime::mugration::pipeline::{self, MugrationInput, MugrationOutput, MugrationParams};
use treetime::progress::RunWarningKind;
use treetime::progress::{LogSink, StageSink};
use treetime::{progress_info, progress_warn};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pair_by_name::pair_by_name;
use treetime_io::discrete_states_csv::discrete_attrs_read_file;
use treetime_utils::io::json::{JsonPretty, json_write_file};

const UNMATCHED_NAMES_SHOWN: usize = 10;

pub fn run_mugration(
  mugration_args: &TreetimeMugrationArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  cancel.check()?;
  stages.report("Reading input", 0.0, "");
  let parse = read_input_tree(&mugration_args.tree, log)?;
  let names = parse.names();
  let graph: Graph = parse.graph;
  let branch_lengths = parse.branch_lengths;

  let resolved = mugration_args.resolve_outputs()?;

  let attribute_column = Some(mugration_args.attribute());

  let (rows, _attr_name) = discrete_attrs_read_file::<String>(
    mugration_args.metadata(),
    &mugration_args.metadata_id.metadata_delimiters,
    &mugration_args.metadata_id.metadata_id_columns,
    None,
    attribute_column,
    |s| Ok(s.to_owned()),
  )?;
  let MugrationTraits {
    traits,
    observed_values,
  } = pair_traits(rows, &graph, &names, mugration_args.metadata(), log);

  let weights = if let Some(weights_filepath) = &mugration_args.weights {
    let (map, _) = discrete_attrs_read_file::<f64>(
      weights_filepath,
      &mugration_args.metadata_id.metadata_delimiters,
      &[],
      attribute_column,
      Some("weight"),
      |s| Ok(s.parse::<f64>()?),
    )?;
    Some(map.into_iter().collect::<BTreeMap<String, f64>>())
  } else {
    None
  };

  cancel.check()?;
  stages.report("Mugration inference", 0.3, "");
  let params = MugrationParams {
    missing_data: mugration_args.missing_data.clone(),
    pc: mugration_args.pc,
    missing_weights_threshold: mugration_args.missing_weights_threshold,
    iterations: mugration_args.iterations,
    sampling_bias_correction: mugration_args.sampling_bias_correction,
    smooth_initial_pi: mugration_args.smooth_initial_pi,
    filter_uninformative_root: mugration_args.filter_uninformative_root,
  };
  let input = MugrationInput {
    graph,
    traits,
    observed_values,
    weights,
    branch_lengths: branch_lengths.clone(),
  };
  let mut output = pipeline::run(&params, input, &names, cancel, log).map_err(|err| err.into_report())?;

  let topology_order = mugration_args
    .topology_order
    .resolve_topology_order(&output.graph, &names, None)?;
  topology_order.apply(&mut output.graph, &names, &branch_lengths)?;
  stages.report("Writing output", 0.8, "");

  if let Some(path) = resolved.path(OutputSelection::Gtr) {
    let gtr_output = GtrOutput::builder()
      .gtr(&output.gtr)
      .model_name(GtrModelName::Infer)
      .attribute(mugration_args.attribute())
      .states(output.states.iter().map(ToOwned::to_owned).collect())
      .build();
    json_write_file(path, &gtr_output, JsonPretty(true))?;
  }

  write_mugration_trees(
    &output,
    &names,
    &branch_lengths,
    mugration_args.attribute(),
    &resolved,
    log,
  )?;

  stages.report("Done", 1.0, "");
  Ok(())
}

pub(crate) fn pair_traits(
  rows: Vec<(String, String)>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  path: &Path,
  log: &dyn LogSink,
) -> MugrationTraits {
  let pairing = pair_by_name(graph.get_leaves().map(|leaf| leaf.key()), names, rows);
  warn_duplicate_names(
    log,
    RunWarningKind::DuplicateMetadataNames,
    &format!("The metadata '{}' has more than one row named", path.display()),
    "TreeTime uses the first row of each name.",
    &pairing.duplicate_entry_names,
  );
  if !pairing.unmatched.is_empty() {
    let shown = pairing
      .unmatched
      .iter()
      .take(UNMATCHED_NAMES_SHOWN)
      .map(|(name, _)| name)
      .join(", ");
    let suffix = if pairing.unmatched.len() > UNMATCHED_NAMES_SHOWN {
      "..."
    } else {
      ""
    };
    progress_warn!(
      log,
      "Mugration: {} metadata names not present in tree: {shown}{suffix}",
      pairing.unmatched.len()
    );
  }
  let observed_values = pairing
    .by_node
    .values()
    .chain(pairing.unmatched.iter().map(|(_, value)| value))
    .cloned()
    .collect();
  MugrationTraits {
    traits: pairing.by_node,
    observed_values,
  }
}

pub(crate) struct MugrationTraits {
  pub(crate) traits: BTreeMap<GraphNodeKey, String>,
  pub(crate) observed_values: BTreeSet<String>,
}

fn write_mugration_trees(
  output: &MugrationOutput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  attribute: &str,
  resolved: &ResolvedOutputs,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let annotated = AnnotatedGraph {
    graph: &output.graph,
    names,
    divergence_branch_lengths: branch_lengths,
    time_branch_lengths: None,
    divergence: Divergence::CumulativeBranchLength,
    sequences: None,
    dates: None,
    traits: Some(TreeTraits {
      attribute,
      states: &output.states,
      values: &output.reconstructed_traits,
      profiles: &output.confidences,
    }),
  };
  write_graph_outputs(&annotated, &resolved.tree_outputs)?;
  let Some(tree) = tree_view_for_outputs(&annotated, resolved)? else {
    return Ok(());
  };
  write_tree_outputs(&tree, &resolved.tree_outputs, CommandKind::Mugration, log)?;
  if let Some(path) = resolved.path(OutputSelection::TraitsCsv) {
    write_traits_csv(&tree, path)?;
  }
  if let Some(path) = resolved.path(OutputSelection::ConfidenceCsv) {
    write_trait_confidence_csv(&tree, path)?;
  }
  if let Some(path) = resolved.path(OutputSelection::AugurNodeData) {
    write_augur_node_data_traits(&tree, &output.gtr, path)?;
    progress_info!(log, "Wrote augur node data JSON to {}", path.display());
  }
  Ok(())
}
