use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::commands::mugration::augur_node_data::write_augur_node_data_json;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use app_output::discrete_trait_comment::DiscreteTraitCommentProvider;
use app_output::mugration_result::MugrationResult;
use app_output::mugration_tree_output::write_mugration_tree_outputs;
use app_output::output_plan::OutputSelection;
use eyre::Report;
use log::info;
use std::collections::BTreeMap;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput, write_gtr_json};
use treetime::make_report;
use treetime::mugration::pipeline::{self, MugrationInput, MugrationParams};
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
  let mut output = pipeline::run(&params, input, &names, cancel).map_err(|err| err.into_report())?;

  let topology_order = mugration_args
    .topology_order
    .resolve_topology_order(&output.graph, &names, None)?;
  topology_order.apply(&mut output.graph, &names, &branch_lengths)?;
  progress.report("Writing output", 0.8, "");

  // Project the canonical core output into the serializable per-file views the writers consume. The
  // reconstructed value maps and the tree topology stay in `output`; the writers read the maps from it
  // and the per-node/per-edge metadata from `result`.
  let result = MugrationResult::new(
    &output,
    &confidences,
    &names,
    &branch_lengths,
    mugration_args.attribute(),
  );

  if !resolved.tree_outputs.is_empty() {
    let provider = DiscreteTraitCommentProvider::new(&output.reconstructed_traits, mugration_args.attribute());
    let providers = CommentProviders::new().with(&provider);
    write_mugration_tree_outputs(
      &output,
      &result.nodes,
      &branch_lengths,
      mugration_args.attribute(),
      &resolved.tree_outputs,
      &providers,
    )?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    let gtr_output = GtrOutput::new(&output.gtr, GtrModelName::Infer)
      .with_discrete_states(mugration_args.attribute(), output.states.iter());
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
    write_augur_node_data_json(&result, &output, path)?;
    info!("Wrote augur node data JSON to {}", path.display());
  }

  progress.report("Done", 1.0, "");
  Ok(result)
}
