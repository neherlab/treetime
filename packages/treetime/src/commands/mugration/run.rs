use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::commands::mugration::comment_provider::PartitionCommentProvider;
use crate::mugration::input::MugrationInput;
use crate::mugration::mugration::execute_mugration;
use crate::mugration::result::{MugrationGtrOutput, MugrationResult, MugrationTraitsOutput};
use crate::partition::marginal_discrete::PartitionMarginalDiscrete;
use crate::payload::ancestral::GraphAncestral;
use crate::make_report;
use eyre::Report;
use log::info;
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use treetime_io::discrete_states_csv::read_discrete_attrs;
use treetime_io::nex::{NexWriteOptions, nex_write_str_with};
use treetime_io::nwk::CommentProviders;
use treetime_io::nwk::nwk_read_file;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn parse_mugration_input(args: &TreetimeMugrationArgs) -> Result<MugrationInput, Report> {
  let TreetimeMugrationArgs {
    tree,
    attribute,
    states,
    weights,
    name_column,
    pc,
    missing_data,
    missing_weights_threshold,
    ..
  } = args;

  let tree_path = tree.as_ref().ok_or_else(|| make_report!("Tree file is required"))?;
  let graph: GraphAncestral = nwk_read_file(tree_path)?;

  let (attr_values, _attr_name) =
    read_discrete_attrs::<String>(states, name_column, &Some(attribute.clone()), |s| Ok(s.to_owned()))?;
  let traits: BTreeMap<String, String> = attr_values.into_iter().collect();

  let weights_map = if let Some(weights_filepath) = weights {
    let (map, _) = read_discrete_attrs::<f64>(
      weights_filepath,
      &Some(attribute.clone()),
      &Some("weight".to_owned()),
      |s| Ok(s.parse::<f64>()?),
    )?;
    Some(map.into_iter().collect())
  } else {
    None
  };

  Ok(MugrationInput {
    graph,
    traits,
    attribute: attribute.clone(),
    weights: weights_map,
    missing_data: missing_data.clone(),
    pc: *pc,
    missing_weights_threshold: *missing_weights_threshold,
    iterations: args.iterations,
    sampling_bias_correction: args.sampling_bias_correction,
  })
}

pub fn run_mugration(mugration_args: &TreetimeMugrationArgs) -> Result<(), Report> {
  let outdir = &mugration_args.outdir;
  fs::create_dir_all(outdir)?;

  let input = parse_mugration_input(mugration_args)?;
  let result = execute_mugration(input)?;

  write_annotated_tree(&result.graph, &result.partition, &result.traits, outdir)?;
  write_gtr_json_file(&result.gtr, outdir)?;

  if let Some(confidence_path) = &mugration_args.confidence {
    write_confidence_csv(&result, confidence_path)?;
  }

  info!("Mugration: wrote output to {}", outdir.display());
  Ok(())
}

fn write_annotated_tree(
  graph: &GraphAncestral,
  partition: &PartitionMarginalDiscrete,
  traits: &MugrationTraitsOutput,
  outdir: &Path,
) -> Result<(), Report> {
  let provider = PartitionCommentProvider::new(partition, &traits.attribute);
  let providers = CommentProviders::new().with(&provider);
  let nexus = nex_write_str_with(graph, &NexWriteOptions::default(), &providers)?;
  fs::write(outdir.join("annotated_tree.nexus"), format!("{nexus}\n"))?;

  fs::write(outdir.join("traits.csv"), traits.render_csv())?;

  Ok(())
}

fn write_gtr_json_file(gtr: &MugrationGtrOutput, outdir: &Path) -> Result<(), Report> {
  json_write_file(outdir.join("gtr.json"), gtr, JsonPretty(true))
}

fn write_confidence_csv(result: &MugrationResult, output_path: &Path) -> Result<(), Report> {
  fs::write(output_path, result.confidence.render_csv())?;
  Ok(())
}
