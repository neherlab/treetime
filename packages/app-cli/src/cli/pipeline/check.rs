use crate::cli::pipeline::inputs::labeled_input_paths;
use crate::cli::pipeline::resolve::{ResolvedPipeline, ResolvedStep};
use crate::cli::pipeline::runner::select_steps;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Write};

pub fn print_pipeline_plan(pipeline: &ResolvedPipeline, selected: Option<&BTreeSet<String>>) -> Result<(), Report> {
  let steps = select_steps(pipeline, selected)?;
  let producers = producers_by_path(pipeline);
  let mut out = io::stdout().lock();

  if let Some(workdir) = &pipeline.workdir {
    writeln!(out, "workdir: {}", workdir.display())?;
  }
  writeln!(out, "steps ({}):", steps.len())?;

  for step in steps {
    writeln!(out, "  - {} ({})", step.name, step.command.tag())?;
    print_inputs(&mut out, step, &producers)?;
    print_outputs(&mut out, step)?;
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

fn print_inputs(out: &mut impl Write, step: &ResolvedStep, producers: &BTreeMap<String, String>) -> Result<(), Report> {
  for (label, path) in labeled_input_paths(&step.command)? {
    match producers.get(&path) {
      Some(producer) if producer != &step.name => writeln!(out, "    {label}: {path} (from step {producer})")?,
      _ => writeln!(out, "    {label}: {path}")?,
    }
  }
  Ok(())
}

fn print_outputs(out: &mut impl Write, step: &ResolvedStep) -> Result<(), Report> {
  if let Some(dir) = &step.outputs.output_all {
    writeln!(out, "    output dir: {}", dir.display())?;
  }
  let produced: Vec<String> = step
    .outputs
    .by_selection
    .values()
    .flatten()
    .map(|path| path.to_string_lossy().into_owned())
    .collect();
  if !produced.is_empty() {
    writeln!(out, "    produces:")?;
    for path in produced {
      writeln!(out, "      {path}")?;
    }
  }
  Ok(())
}
