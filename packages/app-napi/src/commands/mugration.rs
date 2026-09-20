use crate::commands::support::{default_topology_order, output_plan};
use app_output::augur_node_data_mugration::write_augur_node_data_json;
use app_output::discrete_trait_comment::DiscreteTraitCommentProvider;
use app_output::mugration_result::MugrationResult;
use app_output::mugration_tree_output::write_mugration_tree_outputs;
use app_output::output_plan::{CommandKind, OutputSelection};
use eyre::Report;
use log::info;
use serde::Deserialize;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime::cancel::Cancel;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput, write_gtr_json};
use treetime::make_report;
use treetime::mugration::pipeline::{self, MugrationInput, MugrationParams};
use treetime::progress::ProgressSink;
use treetime_io::csv::{default_metadata_delimiters, default_name_candidates};
use treetime_io::discrete_states_csv::read_discrete_attrs;
use treetime_io::nwk::{CommentProviders, nwk_read_file};
use treetime_utils::io::file::create_file_or_stdout;

#[derive(Debug, SmartDefault, Deserialize)]
#[serde(default)]
pub struct MugrationArgs {
  pub tree: Option<String>,
  #[default(_code = r#""country".to_owned()"#)]
  pub attribute: String,
  pub states: String,
  pub weights: Option<String>,
  pub name_column: Option<String>,
  pub confidence: Option<String>,
  pub pc: Option<f64>,
  #[default(_code = r#""?".to_owned()"#)]
  pub missing_data: String,
  #[default = 0.5]
  pub missing_weights_threshold: f64,
  #[default = 5]
  pub iterations: usize,
  pub sampling_bias_correction: Option<f64>,
  pub outdir: String,
}

pub fn run_mugration(
  args: &MugrationArgs,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<MugrationResult, Report> {
  cancel.check()?;
  progress.report("Reading input", 0.0, "");
  let tree_path = args
    .tree
    .as_ref()
    .ok_or_else(|| make_report!("Tree file is required"))?;
  let parse = nwk_read_file(Path::new(tree_path))?;
  let confidences = parse.confidences();
  let names = parse.names();
  let graph = parse.graph;
  let branch_lengths = parse.branch_lengths;

  let non_tree_overrides: BTreeMap<OutputSelection, PathBuf> = args
    .confidence
    .as_ref()
    .map(|path| BTreeMap::from([(OutputSelection::ConfidenceCsv, PathBuf::from(path))]))
    .unwrap_or_default();
  let resolved = output_plan(CommandKind::Mugration, Path::new(&args.outdir), non_tree_overrides)?;

  let attribute_column = Some(args.attribute.clone());
  let delimiters = default_metadata_delimiters();
  let id_columns = args
    .name_column
    .clone()
    .map_or_else(default_name_candidates, |col| vec![col]);

  let (attr_values, _attr_name) = read_discrete_attrs::<String>(
    Path::new(&args.states),
    &delimiters,
    &id_columns,
    &None,
    &attribute_column,
    |s| Ok(s.to_owned()),
  )?;
  let traits: BTreeMap<String, String> = attr_values.into_iter().collect();

  let weights = if let Some(weights_filepath) = &args.weights {
    let (map, _) = read_discrete_attrs::<f64>(
      Path::new(weights_filepath),
      &delimiters,
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
    missing_data: args.missing_data.clone(),
    pc: args.pc,
    missing_weights_threshold: args.missing_weights_threshold,
    iterations: args.iterations,
    sampling_bias_correction: args.sampling_bias_correction,
    smooth_initial_pi: false,
    filter_uninformative_root: false,
  };
  let input = MugrationInput {
    graph,
    traits,
    weights,
    branch_lengths: branch_lengths.clone(),
  };
  let mut output = pipeline::run(&params, input, &names, cancel).map_err(|err| err.into_report())?;

  default_topology_order().apply(&mut output.graph, &names, &branch_lengths)?;
  progress.report("Writing output", 0.8, "");

  let result = MugrationResult::new(&output, &confidences, &names, &branch_lengths, &args.attribute);

  if !resolved.tree_outputs.is_empty() {
    let provider = DiscreteTraitCommentProvider::new(&output.reconstructed_traits, &args.attribute);
    let providers = CommentProviders::new().with(&provider);
    write_mugration_tree_outputs(
      &output,
      &result.nodes,
      &branch_lengths,
      &args.attribute,
      &resolved.tree_outputs,
      &providers,
    )?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) {
    let gtr_output =
      GtrOutput::new(&output.gtr, GtrModelName::Infer).with_discrete_states(&args.attribute, output.states.iter());
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
