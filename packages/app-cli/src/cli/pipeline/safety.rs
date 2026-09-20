use crate::cli::pipeline::inputs::input_paths;
use crate::cli::pipeline::resolve::{ResolvedPipeline, ResolvedStep, is_template_path};
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use treetime_utils::make_error;

const STDIO_PATH: &str = "-";

pub fn validate_plan(pipeline: &ResolvedPipeline, selected: Option<&BTreeSet<String>>) -> Result<(), Report> {
  reject_stdio_paths(pipeline)?;
  reject_output_collisions(pipeline)?;
  reject_self_truncation(pipeline)?;
  reject_missing_upstream(pipeline, selected)?;
  Ok(())
}

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

  #[test]
  fn test_safety_output_collision_rejected() {
    let pipeline = resolve(json!({
      "steps": [
        { "name": "a", "ancestral": { "tree": "in.nwk", "output_tree_nwk": "shared.nwk" } },
        { "name": "b", "ancestral": { "tree": "in.nwk", "output_tree_nwk": "shared.nwk" } }
      ]
    }));
    let result = validate_plan(&pipeline, None);
    assert_error!(
      result,
      "steps `a` and `b` both write `shared.nwk`; each output path must be written by only one step"
    );
  }

  #[test]
  fn test_safety_self_truncation_rejected() {
    let pipeline = resolve(json!({
      "steps": [
        { "name": "a", "ancestral": { "tree": "same.nwk", "output_tree_nwk": "same.nwk" } }
      ]
    }));
    let result = validate_plan(&pipeline, None);
    assert_error!(
      result,
      "step `a` reads and writes the same file `same.nwk`, which would truncate the input"
    );
  }

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
