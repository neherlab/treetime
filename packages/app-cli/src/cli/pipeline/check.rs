use crate::cli::pipeline::inputs::labeled_input_paths;
use crate::cli::pipeline::resolve::{ResolvedPipeline, ResolvedStep};
use crate::cli::pipeline::runner::select_steps;
use eyre::{Report, WrapErr};
use itertools::izip;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Write};
use std::path::Path;

pub(crate) fn print_pipeline_plan(
  pipeline: &ResolvedPipeline,
  selected: Option<&BTreeSet<String>>,
) -> Result<(), Report> {
  let steps = select_steps(pipeline, selected)?;
  let producers = producers_by_path(pipeline);
  let step_inputs = steps
    .iter()
    .map(|step| labeled_input_paths(&step.command))
    .collect::<Result<Vec<_>, Report>>()?;

  write_pipeline_plan(
    &mut io::stdout().lock(),
    pipeline.workdir.as_deref(),
    &steps,
    &step_inputs,
    &producers,
  )
  .wrap_err("When writing the pipeline plan to standard output")
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

fn write_pipeline_plan(
  out: &mut impl Write,
  workdir: Option<&Path>,
  steps: &[&ResolvedStep],
  step_inputs: &[Vec<(&'static str, String)>],
  producers: &BTreeMap<String, String>,
) -> io::Result<()> {
  if let Some(workdir) = workdir {
    writeln!(out, "workdir: {}", workdir.display())?;
  }
  writeln!(out, "steps ({}):", steps.len())?;

  for (step, inputs) in izip!(steps, step_inputs) {
    writeln!(out, "  - {} ({})", step.name, step.command.tag())?;
    write_inputs(out, step, inputs, producers)?;
    write_outputs(out, step)?;
  }
  Ok(())
}

fn write_inputs(
  out: &mut impl Write,
  step: &ResolvedStep,
  inputs: &[(&'static str, String)],
  producers: &BTreeMap<String, String>,
) -> io::Result<()> {
  for (label, path) in inputs {
    match producers.get(path) {
      Some(producer) if producer != &step.name => writeln!(out, "    {label}: {path} (from step {producer})")?,
      _ => writeln!(out, "    {label}: {path}")?,
    }
  }
  Ok(())
}

fn write_outputs(out: &mut impl Write, step: &ResolvedStep) -> io::Result<()> {
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
