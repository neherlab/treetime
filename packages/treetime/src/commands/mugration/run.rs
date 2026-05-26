use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::gtr::get_gtr::{GtrModelName, GtrOutput};
use crate::make_report;
use crate::mugration::mugration::execute_mugration;
use crate::mugration::result::MugrationResult;
use crate::partition::marginal_discrete::DiscreteCommentProvider;
use crate::partition::traits::HasGtr;
use crate::payload::ancestral::GraphAncestral;
use eyre::Report;
use log::info;
use std::collections::BTreeMap;
use std::fs;
use treetime_io::discrete_states_csv::read_discrete_attrs;
use treetime_io::graph::write_graph_files_with;
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn run_mugration(
  mugration_args: &TreetimeMugrationArgs,
  progress: &dyn crate::progress::ProgressSink,
) -> Result<MugrationResult, Report> {
  let outdir = &mugration_args.outdir;
  fs::create_dir_all(outdir)?;

  progress.check_cancelled()?;
  progress.report("Reading input", 0.0, "");
  let tree_path = mugration_args
    .tree
    .as_ref()
    .ok_or_else(|| make_report!("Tree file is required"))?;
  let graph: GraphAncestral = nwk_read_file(tree_path)?;

  let (attr_values, _attr_name) = read_discrete_attrs::<String>(
    &mugration_args.states,
    &mugration_args.name_column,
    &Some(mugration_args.attribute.clone()),
    |s| Ok(s.to_owned()),
  )?;
  let traits: BTreeMap<String, String> = attr_values.into_iter().collect();

  let weights = if let Some(weights_filepath) = &mugration_args.weights {
    let (map, _) = read_discrete_attrs::<f64>(
      weights_filepath,
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
  let result = execute_mugration(
    graph,
    &traits,
    &mugration_args.attribute,
    weights.as_ref(),
    &mugration_args.missing_data,
    mugration_args.pc,
    mugration_args.missing_weights_threshold,
    mugration_args.iterations,
    mugration_args.sampling_bias_correction,
  )?;

  progress.report("Writing output", 0.8, "");
  let provider = DiscreteCommentProvider::new(&result.partition, &result.traits.attribute);
  let providers = CommentProviders::new().with(&provider);
  write_graph_files_with(outdir, "annotated_tree", &result.graph, &providers)?;

  let gtr_output = GtrOutput::new(result.partition.gtr(), GtrModelName::Infer)
    .with_discrete_states(&result.traits.attribute, result.partition.states.iter());
  json_write_file(outdir.join("gtr.json"), &gtr_output, JsonPretty(true))?;

  fs::write(outdir.join("traits.csv"), result.traits.render_csv())?;

  if let Some(confidence_path) = &mugration_args.confidence {
    fs::write(confidence_path, result.confidence.render_csv())?;
  }

  progress.report("Done", 1.0, "");
  info!("Mugration: wrote output to {}", outdir.display());
  Ok(result)
}
