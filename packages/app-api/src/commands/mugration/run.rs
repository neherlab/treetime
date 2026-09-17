use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::commands::mugration::augur_node_data::write_augur_node_data_json;
use crate::commands::shared::output::OutputSelection;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use app_output::mugration_tree_output::write_mugration_tree_outputs;
use eyre::Report;
use log::info;
use std::collections::BTreeMap;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput, write_gtr_json};
use treetime::make_report;
use treetime::mugration::mugration::execute_mugration;
use treetime::mugration::result::MugrationResult;
use treetime::partition::marginal::discrete::comment::DiscreteTraitCommentProvider;
use treetime_graph::graph::Graph;
use treetime_io::discrete_states_csv::read_discrete_attrs;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_utils::io::file::create_file_or_stdout;

pub fn run_mugration(
  mugration_args: &TreetimeMugrationArgs,
  cancel: &dyn treetime::cancel::Cancel,
  progress: &dyn treetime::progress::ProgressSink,
) -> Result<MugrationResult, Report> {
  cancel.check()?;
  progress.report("Reading input", 0.0, "");
  let tree_path = mugration_args
    .tree
    .as_ref()
    .ok_or_else(|| make_report!("Tree file is required"))?;
  let parse = nwk_read_file(tree_path)?;
  let confidences = parse.confidences();
  let names = parse.names();
  let graph: Graph = parse.graph;
  let branch_lengths = parse.branch_lengths;

  let resolved = mugration_args.resolve_outputs()?;

  // The attribute names the metadata column to read discrete states from.
  let attribute_column = Some(mugration_args.attribute().to_owned());

  let (attr_values, _attr_name) = read_discrete_attrs::<String>(
    mugration_args.metadata(),
    &mugration_args.metadata_id.metadata_delimiters,
    &mugration_args.metadata_id.metadata_id_columns,
    &None,
    &attribute_column,
    |s| Ok(s.to_owned()),
  )?;
  let traits: BTreeMap<String, String> = attr_values.into_iter().collect();

  let weights = if let Some(weights_filepath) = &mugration_args.weights {
    let (map, _) = read_discrete_attrs::<f64>(
      weights_filepath,
      &mugration_args.metadata_id.metadata_delimiters,
      &[],
      &attribute_column,
      &Some("weight".to_owned()),
      |s| Ok(s.parse::<f64>()?),
    )?;
    Some(map.into_iter().collect::<BTreeMap<String, f64>>())
  } else {
    None
  };

  cancel.check()?;
  progress.report("Mugration inference", 0.3, "");
  // The output value maps are gathered off the pipeline-local partition inside `execute_mugration`,
  // before the partition-less graph data slot is built; the tree, Newick-comment, augur, and GTR
  // writers read the maps here.
  let (mut result, maps) = execute_mugration(
    graph,
    &confidences,
    &names,
    &branch_lengths,
    &traits,
    mugration_args.attribute(),
    weights.as_ref(),
    &mugration_args.missing_data,
    mugration_args.pc,
    mugration_args.missing_weights_threshold,
    mugration_args.iterations,
    mugration_args.sampling_bias_correction,
    mugration_args.smooth_initial_pi,
    mugration_args.filter_uninformative_root,
    cancel,
  )?;

  let topology_order = mugration_args
    .topology_order
    .resolve_topology_order(&result.graph, &names, None)?;
  topology_order.apply(&mut result.graph, &names, &branch_lengths)?;
  progress.report("Writing output", 0.8, "");

  if !resolved.tree_outputs.is_empty() {
    let provider = DiscreteTraitCommentProvider::new(&maps.reconstructed_traits, &result.traits.attribute);
    let providers = CommentProviders::new().with(&provider);
    write_mugration_tree_outputs(
      &result.graph,
      &result.nodes,
      &branch_lengths,
      &maps,
      &result.traits.attribute,
      &resolved.tree_outputs,
      &providers,
    )?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    let gtr_output =
      GtrOutput::new(&maps.gtr, GtrModelName::Infer).with_discrete_states(&result.traits.attribute, maps.states.iter());
    write_gtr_json(&gtr_output, path)?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::TraitsCsv) {
    let mut f = create_file_or_stdout(path)?;
    std::io::Write::write_all(&mut f, result.traits.render_csv().as_bytes())?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ConfidenceCsv) {
    let mut f = create_file_or_stdout(path)?;
    std::io::Write::write_all(&mut f, result.confidence.render_csv().as_bytes())?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    write_augur_node_data_json(&result, &maps, path)?;
    info!("Wrote augur node data JSON to {}", path.display());
  }

  progress.report("Done", 1.0, "");
  Ok(result)
}
