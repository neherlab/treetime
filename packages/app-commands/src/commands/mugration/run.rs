use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeTraits};
use app_output::augur_node_data_traits::write_augur_node_data_traits;
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
use app_output::trait_tables::{write_trait_confidence_csv, write_traits_csv};
use app_output::tree_output::{tree_view_for_outputs, write_graph_outputs, write_tree_outputs};
use eyre::Report;
use std::collections::BTreeMap;
use treetime::cancel::Cancel;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput};
use treetime::mugration::pipeline::{self, MugrationInput, MugrationOutput, MugrationParams};
use treetime::progress::{LogSink, StageSink};
use treetime::progress_info;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::discrete_states_csv::discrete_attrs_read_file;
use treetime_io::nwk::nwk_read_file;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn run_mugration(
  mugration_args: &TreetimeMugrationArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  cancel.check()?;
  stages.report("Reading input", 0.0, "");
  let parse = nwk_read_file(&mugration_args.tree)?;
  let confidences = parse.confidences();
  let names = parse.names();
  let graph: Graph = parse.graph;
  let branch_lengths = parse.branch_lengths;

  let resolved = mugration_args.resolve_outputs()?;

  let attribute_column = Some(mugration_args.attribute());

  let (attr_values, _attr_name) = discrete_attrs_read_file::<String>(
    mugration_args.metadata(),
    &mugration_args.metadata_id.metadata_delimiters,
    &mugration_args.metadata_id.metadata_id_columns,
    None,
    attribute_column,
    |s| Ok(s.to_owned()),
  )?;
  let traits: BTreeMap<String, String> = attr_values.into_iter().collect();

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
    weights,
    branch_lengths: branch_lengths.clone(),
  };
  let mut output = pipeline::run(&params, input, &names, cancel, log).map_err(|err| err.into_report())?;

  let topology_order = mugration_args
    .topology_order
    .resolve_topology_order(&output.graph, &names, None)?;
  topology_order.apply(&mut output.graph, &names, &branch_lengths)?;
  stages.report("Writing output", 0.8, "");

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
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
    &confidences,
    mugration_args.attribute(),
    &resolved,
    log,
  )?;

  stages.report("Done", 1.0, "");
  Ok(())
}

fn write_mugration_trees(
  output: &MugrationOutput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  branch_support: &BTreeMap<GraphNodeKey, Option<f64>>,
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
    branch_support: Some(branch_support),
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
  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::TraitsCsv) {
    write_traits_csv(&tree, path)?;
  }
  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ConfidenceCsv) {
    write_trait_confidence_csv(&tree, path)?;
  }
  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    write_augur_node_data_traits(&tree, &output.gtr, path)?;
    progress_info!(log, "Wrote augur node data JSON to {}", path.display());
  }
  Ok(())
}
