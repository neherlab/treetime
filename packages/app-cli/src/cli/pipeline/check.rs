use crate::cli::pipeline::inputs::labeled_input_paths;
use crate::cli::pipeline::resolve::{ResolvedPipeline, ResolvedStep};
use crate::cli::pipeline::runner::select_steps;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};

pub fn print_pipeline_plan(pipeline: &ResolvedPipeline, selected: Option<&BTreeSet<String>>) -> Result<(), Report> {
  let steps = select_steps(pipeline, selected)?;
  let producers = producers_by_path(pipeline);

  if let Some(workdir) = &pipeline.workdir {
    println!("workdir: {}", workdir.display());
  }
  println!("steps ({}):", steps.len());

  for step in steps {
    println!("  - {} ({})", step.name, step.command.tag());
    print_inputs(step, &producers)?;
    print_outputs(step);
  }
  Ok(())
}

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

fn print_inputs(step: &ResolvedStep, producers: &BTreeMap<String, String>) -> Result<(), Report> {
  for (label, path) in labeled_input_paths(&step.command)? {
    match producers.get(&path) {
      Some(producer) if producer != &step.name => println!("    {label}: {path} (from step {producer})"),
      _ => println!("    {label}: {path}"),
    }
  }
  Ok(())
}

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
