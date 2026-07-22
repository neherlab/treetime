use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::commands::mugration::augur_node_data::write_augur_node_data_json;
use crate::commands::shared::output::{CommandKind, OutputSelection};
use crate::commands::shared::tree_output::TreeOutputAdapter;
use crate::gtr::get_gtr::{GtrModelName, GtrOutput, write_gtr_json};
use crate::make_report;
use crate::mugration::mugration::execute_mugration;
use crate::mugration::result::MugrationResult;
use crate::partition::marginal_discrete::DiscreteCommentProvider;
use crate::partition::traits::HasGtr;
use crate::payload::ancestral::GraphAncestral;
use eyre::Report;
use log::info;
use std::collections::BTreeMap;
use treetime_io::discrete_states_csv::read_discrete_attrs;
use treetime_io::graph::write_tree_outputs;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_utils::io::file::create_file_or_stdout;

pub fn run_mugration(
  mugration_args: &TreetimeMugrationArgs,
  progress: &dyn crate::progress::ProgressSink,
) -> Result<MugrationResult, Report> {
  progress.check_cancelled()?;
  progress.report("Reading input", 0.0, "");
  let tree_path = mugration_args
    .tree
    .as_ref()
    .ok_or_else(|| make_report!("Tree file is required"))?;
  let graph: GraphAncestral = nwk_read_file(tree_path)?;

  let selection: Vec<OutputSelection> = mugration_args
    .output_selection
    .iter()
    .copied()
    .map(OutputSelection::from)
    .collect();
  let resolved = mugration_args.output.resolve(
    CommandKind::Mugration,
    &selection,
    &[
      (
        OutputSelection::AugurNodeData,
        mugration_args.output_augur_node_data.as_deref(),
      ),
      (OutputSelection::Gtr, mugration_args.output_gtr.as_deref()),
      (
        OutputSelection::ConfidenceCsv,
        mugration_args.output_confidence_csv.as_deref(),
      ),
      (OutputSelection::TraitsCsv, mugration_args.output_traits_csv.as_deref()),
    ],
  )?;

  let (attr_values, _attr_name) = read_discrete_attrs::<String>(
    &mugration_args.metadata,
    &mugration_args.metadata_id.metadata_id_columns,
    &None,
    &Some(mugration_args.attribute.clone()),
    |s| Ok(s.to_owned()),
  )?;
  let traits: BTreeMap<String, String> = attr_values.into_iter().collect();

  let weights = if let Some(weights_filepath) = &mugration_args.weights {
    let (map, _) = read_discrete_attrs::<f64>(
      weights_filepath,
      &[],
      &Some(mugration_args.attribute.clone()),
      &Some("weight".to_owned()),
      |s| Ok(s.parse::<f64>()?),
    )?;
    Some(map.into_iter().collect::<BTreeMap<String, f64>>())
  } else {
    None
  };

  progress.check_cancelled()?;
  progress.report("Mugration inference", 0.3, "");
  let mut result = execute_mugration(
    graph,
    &traits,
    &mugration_args.attribute,
    weights.as_ref(),
    &mugration_args.missing_data,
    mugration_args.pc,
    mugration_args.missing_weights_threshold,
    mugration_args.iterations,
    mugration_args.sampling_bias_correction,
    mugration_args.smooth_initial_pi,
    mugration_args.filter_uninformative_root,
  )?;

  let topology_order = mugration_args
    .topology_order
    .resolve_topology_order(&result.graph, None)?;
  topology_order.apply(&mut result.graph)?;
  resolved.prepare()?;

  progress.report("Writing output", 0.8, "");

  if !resolved.tree_outputs.is_empty() {
    let provider = DiscreteCommentProvider::new(&result.graph.data().partition, &result.graph.data().traits.attribute);
    let providers = CommentProviders::new().with(&provider);
    write_tree_outputs::<TreeOutputAdapter, _, _, _>(&result.graph, &resolved.tree_outputs, &providers)?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    let gtr_output = GtrOutput::new(result.graph.data().partition.gtr(), GtrModelName::Infer).with_discrete_states(
      &result.graph.data().traits.attribute,
      result.graph.data().partition.states.iter(),
    );
    write_gtr_json(&gtr_output, path)?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::TraitsCsv) {
    let mut f = create_file_or_stdout(path)?;
    std::io::Write::write_all(&mut f, result.graph.data().traits.render_csv().as_bytes())?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ConfidenceCsv) {
    let mut f = create_file_or_stdout(path)?;
    std::io::Write::write_all(&mut f, result.graph.data().confidence.render_csv().as_bytes())?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::AugurNodeData) {
    write_augur_node_data_json(&result, path)?;
    info!("Wrote augur node data JSON to {}", path.display());
  }

  progress.report("Done", 1.0, "");
  Ok(result)
}
