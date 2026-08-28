use crate::cli::pipeline::resolve::{ResolvedPipeline, ResolvedStep};
use crate::cli::pipeline::runner::select_steps;
use eyre::Report;
use serde_json::Value;
use std::collections::BTreeMap;
use std::collections::BTreeSet;

/// Input fields the dry-run plan reports, as (label, path into the args object).
const INPUT_FIELDS: [(&str, &[&str]); 5] = [
  ("tree", &["tree"]),
  ("alignment", &["alignment", "alignment"]),
  ("metadata", &["metadata"]),
  ("weights", &["weights"]),
  ("vcf-reference", &["vcf_reference"]),
];

/// Print the resolved pipeline plan without running anything.
///
/// Shows each step in run order with its resolved input paths (annotated when an input is produced
/// by an earlier step), its output directory, and the files it will produce. Environment values are
/// never printed: a `{{ env.* }}` reference could carry a secret, so only references are shown, and
/// they are already resolved into concrete paths here without exposing the raw variable.
pub fn print_pipeline_plan(pipeline: &ResolvedPipeline, selected: Option<&BTreeSet<String>>) -> Result<(), Report> {
  let steps = select_steps(pipeline, selected)?;
  let producers = producers_by_path(pipeline);

  if let Some(workdir) = &pipeline.workdir {
    println!("workdir: {}", workdir.display());
  }
  println!("steps ({}):", steps.len());

  for step in steps {
    println!("  - {} ({})", step.name, step.command.tag());
    print_inputs(step, &producers);
    print_outputs(step);
  }
  Ok(())
}

/// Map every produced file path back to the step that produces it, for input provenance annotations.
fn producers_by_path(pipeline: &ResolvedPipeline) -> BTreeMap<String, String> {
  let mut producers = BTreeMap::new();
  for step in &pipeline.steps {
    for paths in step.outputs.by_selection.values() {
      for path in paths {
        producers.insert(path.to_string_lossy().into_owned(), step.name.clone());
      }
    }
  }
  producers
}

/// Print a step's resolved input paths, annotating any that an earlier step produces.
fn print_inputs(step: &ResolvedStep, producers: &BTreeMap<String, String>) {
  let args = step.command.args_value();
  for (label, path_keys) in INPUT_FIELDS {
    for path in lookup_paths(&args, path_keys) {
      match producers.get(path) {
        Some(producer) if producer != &step.name => println!("    {label}: {path} (from step {producer})"),
        _ => println!("    {label}: {path}"),
      }
    }
  }
}

/// Print the files a step will produce, grouped in a stable order.
fn print_outputs(step: &ResolvedStep) {
  if let Some(dir) = &step.outputs.output_all {
    println!("    output dir: {}", dir.display());
  }
  let produced: Vec<String> = step
    .outputs
    .by_selection
    .values()
    .flatten()
    .map(|path| path.to_string_lossy().into_owned())
    .collect();
  if !produced.is_empty() {
    println!("    produces:");
    for path in produced {
      println!("      {path}");
    }
  }
}

/// Follow a dotted key path into the args object and return the string paths it holds.
///
/// A leaf may be a single string (`tree`, `metadata`) or a list of strings (`alignment` accepts
/// several files); both are flattened to individual paths.
fn lookup_paths<'a>(args: &'a Value, keys: &[&str]) -> Vec<&'a str> {
  let mut current = args;
  for key in keys {
    let Some(next) = current.get(key) else {
      return Vec::new();
    };
    current = next;
  }
  match current {
    Value::String(path) => vec![path.as_str()],
    Value::Array(items) => items.iter().filter_map(Value::as_str).collect(),
    _ => Vec::new(),
  }
}
