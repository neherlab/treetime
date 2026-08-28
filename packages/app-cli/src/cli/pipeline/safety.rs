use crate::cli::pipeline::inputs::input_paths;
use crate::cli::pipeline::resolve::{ResolvedPipeline, ResolvedStep};
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use treetime_utils::make_error;

/// The stdin/stdout placeholder path, rejected inside a pipeline.
const STDIO_PATH: &str = "-";

/// Placeholders that mark a per-CDS output path template rather than a concrete file.
const CDS_PLACEHOLDERS: [&str; 2] = ["{cds}", "%GENE"];

/// Validate a resolved pipeline before any step runs.
///
/// These checks are cheap and fail the whole run up front, so no partial outputs are written before a
/// structural problem surfaces: two steps writing the same file, a step truncating a file it reads,
/// a selected step reading an absent output of an unselected upstream step, or a `-` (stdin/stdout)
/// path that a single-process multi-step run cannot honor.
pub fn validate_plan(pipeline: &ResolvedPipeline, selected: Option<&BTreeSet<String>>) -> Result<(), Report> {
  reject_stdio_paths(pipeline)?;
  reject_output_collisions(pipeline)?;
  reject_self_truncation(pipeline)?;
  reject_missing_upstream(pipeline, selected)?;
  Ok(())
}

/// Reject `-` used as an input or output path anywhere in the pipeline.
fn reject_stdio_paths(pipeline: &ResolvedPipeline) -> Result<(), Report> {
  for step in &pipeline.steps {
    if input_paths(&step.command).iter().any(|path| path == STDIO_PATH) {
      return make_error!(
        "step `{}` uses `-` (stdin) for an input; a pipeline runs many steps in one process, so steps need real file paths",
        step.name
      );
    }
    if output_paths(step).iter().any(|path| path == STDIO_PATH) {
      return make_error!(
        "step `{}` uses `-` (stdout) for an output; a pipeline runs many steps in one process, so steps need real file paths",
        step.name
      );
    }
  }
  Ok(())
}

/// Reject two steps that resolve to the same output file.
fn reject_output_collisions(pipeline: &ResolvedPipeline) -> Result<(), Report> {
  let mut owner: BTreeMap<String, &str> = BTreeMap::new();
  for step in &pipeline.steps {
    for path in output_paths(step) {
      if let Some(previous) = owner.insert(path.clone(), &step.name) {
        if previous != step.name {
          return make_error!(
            "steps `{previous}` and `{}` both write `{path}`; each output path must be written by only one step",
            step.name
          );
        }
      }
    }
  }
  Ok(())
}

/// Reject a step that reads and writes the same file, which would truncate the input mid-read.
fn reject_self_truncation(pipeline: &ResolvedPipeline) -> Result<(), Report> {
  for step in &pipeline.steps {
    let outputs: BTreeSet<String> = output_paths(step).into_iter().collect();
    for input in input_paths(&step.command) {
      if outputs.contains(&input) {
        return make_error!(
          "step `{}` reads and writes the same file `{input}`, which would truncate the input",
          step.name
        );
      }
    }
  }
  Ok(())
}

/// Reject a selected step that reads an absent output of an unselected upstream step (D8).
///
/// When every step runs there is nothing to check. With a `--steps` subset, a reference to an
/// upstream step's output resolves to that step's path even if the step is skipped; the file must
/// then already exist on disk, otherwise the selected step would read a missing input.
fn reject_missing_upstream(pipeline: &ResolvedPipeline, selected: Option<&BTreeSet<String>>) -> Result<(), Report> {
  let Some(selected) = selected else {
    return Ok(());
  };

  let mut producers: BTreeMap<String, &str> = BTreeMap::new();
  for step in &pipeline.steps {
    for path in output_paths(step) {
      producers.insert(path, &step.name);
    }
  }

  for step in &pipeline.steps {
    if !selected.contains(&step.name) {
      continue;
    }
    for input in input_paths(&step.command) {
      if let Some(producer) = producers.get(&input) {
        if !selected.contains(*producer) && !Path::new(&input).exists() {
          return make_error!(
            "step `{}` reads `{input}`, produced by step `{producer}`, which is not selected and whose output is absent; \
             include `{producer}` in --steps or run it first",
            step.name
          );
        }
      }
    }
  }
  Ok(())
}

/// Concrete output file paths a step writes, excluding per-CDS templates.
fn output_paths(step: &ResolvedStep) -> Vec<String> {
  step
    .outputs
    .by_selection
    .values()
    .flatten()
    .filter(|path| !is_template_path(path))
    .map(|path| path.to_string_lossy().into_owned())
    .collect()
}

/// Whether a path is a per-CDS template rather than a concrete file.
fn is_template_path(path: &Path) -> bool {
  let path = path.to_string_lossy();
  CDS_PLACEHOLDERS.iter().any(|placeholder| path.contains(placeholder))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::cli::pipeline::resolve::{PipelineDoc, resolve_pipeline};
  use maplit::btreeset;
  use serde_json::{Value, json};
  use treetime_utils::assert_error;

  fn resolve(config: Value) -> ResolvedPipeline {
    let doc = PipelineDoc::from_value(config).unwrap();
    resolve_pipeline(&doc, &json!({})).unwrap()
  }

  // Two steps writing the same explicit output file are rejected.
  #[test]
  fn test_safety_output_collision_rejected() {
    let pipeline = resolve(json!({
      "steps": [
        { "name": "a", "ancestral": { "tree": "in.nwk", "output": { "output_tree_nwk": "shared.nwk" } } },
        { "name": "b", "ancestral": { "tree": "in.nwk", "output": { "output_tree_nwk": "shared.nwk" } } }
      ]
    }));
    let result = validate_plan(&pipeline, None);
    assert_error!(
      result,
      "steps `a` and `b` both write `shared.nwk`; each output path must be written by only one step"
    );
  }

  // A step that reads and writes the same file is rejected as a truncation.
  #[test]
  fn test_safety_self_truncation_rejected() {
    let pipeline = resolve(json!({
      "steps": [
        { "name": "a", "ancestral": { "tree": "same.nwk", "output": { "output_tree_nwk": "same.nwk" } } }
      ]
    }));
    let result = validate_plan(&pipeline, None);
    assert_error!(
      result,
      "step `a` reads and writes the same file `same.nwk`, which would truncate the input"
    );
  }

  // A `-` input path is rejected: a multi-step single-process run has one stdin.
  #[test]
  fn test_safety_stdio_input_rejected() {
    let pipeline = resolve(json!({
      "steps": [
        { "name": "a", "ancestral": { "tree": "-", "output_gtr": "g.json" } }
      ]
    }));
    let result = validate_plan(&pipeline, None);
    assert_error!(
      result,
      "step `a` uses `-` (stdin) for an input; a pipeline runs many steps in one process, so steps need real file paths"
    );
  }

  // A selected step referencing an unselected upstream step's absent output is rejected.
  #[test]
  fn test_safety_missing_upstream_rejected() {
    let pipeline = resolve(json!({
      "output_all": "tmp/pipeline-safety-absent",
      "steps": [
        { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv" } },
        { "name": "anc", "ancestral": { "tree": "{{ steps.tt.outputs.nwk }}" } }
      ]
    }));
    let result = validate_plan(&pipeline, Some(&btreeset! { "anc".to_owned() }));
    assert_error!(
      result,
      "step `anc` reads `tmp/pipeline-safety-absent/tt/timetree.nwk`, produced by step `tt`, which is not selected and whose output is absent; include `tt` in --steps or run it first"
    );
  }

  // With every step selected, upstream outputs are produced by the run, so nothing is missing.
  #[test]
  fn test_safety_full_run_has_no_missing_upstream() {
    let pipeline = resolve(json!({
      "output_all": "tmp/pipeline-safety-absent",
      "steps": [
        { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv" } },
        { "name": "anc", "ancestral": { "tree": "{{ steps.tt.outputs.nwk }}" } }
      ]
    }));
    validate_plan(&pipeline, None).unwrap();
  }
}
