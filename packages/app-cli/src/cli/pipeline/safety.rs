use crate::cli::pipeline::inputs::input_paths;
use crate::cli::pipeline::resolve::{ResolvedPipeline, ResolvedStep};
use app_commands::config::cli_rules::cli_rule_diagnostics;
use app_commands::config::source::{ConfigSource, render_and_bail};
use app_output::output_plan::PathPlaceholder;
use eyre::{Report, WrapErr};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use treetime_utils::make_error;

const STDIO_PATH: &str = "-";

pub(crate) fn validate_plan(
  pipeline: &ResolvedPipeline,
  source: &ConfigSource,
  selected: Option<&BTreeSet<String>>,
) -> Result<(), Report> {
  reject_cli_rule_violations(pipeline, source)?;
  reject_stdio_paths(pipeline)?;
  reject_output_collisions(pipeline)?;
  reject_self_truncation(pipeline)?;
  reject_missing_upstream(pipeline, selected)?;
  Ok(())
}

fn reject_cli_rule_violations(pipeline: &ResolvedPipeline, source: &ConfigSource) -> Result<(), Report> {
  for (index, step) in pipeline.steps.iter().enumerate() {
    let command = step.command.command();
    let tag: &str = command.into();
    let diags = cli_rule_diagnostics(command, &step.command.settings()?)?
      .into_iter()
      .map(|mut diag| {
        diag.pointer = diag.pointer.map(|pointer| format!("/steps/{index}/{tag}{pointer}"));
        diag
      })
      .collect();
    render_and_bail(source, "invalid configuration", diags)
      .wrap_err_with(|| format!("in pipeline step `{}`", step.name))?;
  }
  Ok(())
}

fn reject_stdio_paths(pipeline: &ResolvedPipeline) -> Result<(), Report> {
  for step in &pipeline.steps {
    if input_paths(&step.command)?.iter().any(|path| path == STDIO_PATH) {
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
    for input in input_paths(&step.command)? {
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
    for input in input_paths(&step.command)? {
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
    .filter(|path| !PathPlaceholder::Cds.is_in_file_name(path))
    .map(|path| path.to_string_lossy().into_owned())
    .collect()
}

#[cfg(test)]
mod tests {
  const PIPELINE_FOLDER: &str = "/pipeline";

  use super::*;
  use crate::cli::pipeline::resolve::{PipelineDoc, resolve_pipeline};
  use app_commands::config::source::InvalidConfig;
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use treetime_utils::error::report_to_string;
  use treetime_utils::{assert_error, o};

  fn validate(config: Value, selected: Option<&BTreeSet<String>>) -> Result<(), Report> {
    let source = ConfigSource::new("pipeline.json", config.to_string());
    let doc = PipelineDoc::from_value(config).unwrap();
    let pipeline = resolve_pipeline(&doc, &json!({}), Path::new(PIPELINE_FOLDER), None).unwrap();
    validate_plan(&pipeline, &source, selected)
  }

  #[test]
  fn test_safety_output_collision_rejected() {
    let result = validate(
      json!({
        "steps": [
          { "name": "a", "ancestral": { "tree": "in.nwk", "output_tree_nwk": "shared.nwk" } },
          { "name": "b", "ancestral": { "tree": "in.nwk", "output_tree_nwk": "shared.nwk" } }
        ]
      }),
      None,
    );
    assert_error!(
      result,
      "steps `a` and `b` both write `/pipeline/shared.nwk`; each output path must be written by only one step"
    );
  }

  #[test]
  fn test_safety_output_collision_of_a_partition_file_rejected() {
    let result = validate(
      json!({
        "steps": [
          { "name": "a", "ancestral": { "tree": "in.nwk", "alignment": ["x.fasta", "y.fasta"], "output_reconstructed_nuc_fasta": "seq.{partition}.fasta" } },
          { "name": "b", "ancestral": { "tree": "in.nwk", "output_reconstructed_nuc_fasta": "seq.2.fasta" } }
        ]
      }),
      None,
    );
    assert_error!(
      result,
      "steps `a` and `b` both write `/pipeline/seq.2.fasta`; each output path must be written by only one step"
    );
  }

  #[test]
  fn test_safety_self_truncation_rejected() {
    let result = validate(
      json!({
        "steps": [
          { "name": "a", "ancestral": { "tree": "same.nwk", "output_tree_nwk": "same.nwk" } }
        ]
      }),
      None,
    );
    assert_error!(
      result,
      "step `a` reads and writes the same file `/pipeline/same.nwk`, which would truncate the input"
    );
  }

  #[test]
  fn test_safety_stdio_input_rejected() {
    let result = validate(
      json!({
        "steps": [
          { "name": "a", "ancestral": { "tree": "-", "output_gtr": "g.json" } }
        ]
      }),
      None,
    );
    assert_error!(
      result,
      "step `a` uses `-` (stdin) for an input; a pipeline runs many steps in one process, so steps need real file paths"
    );
  }

  #[test]
  fn test_safety_missing_upstream_rejected() {
    let result = validate(
      json!({
        "output_all": "tmp/pipeline-safety-absent",
        "steps": [
          { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv" } },
          { "name": "anc", "ancestral": { "tree": "{{ steps.tt.outputs.nwk }}" } }
        ]
      }),
      Some(&btreeset! { "anc".to_owned() }),
    );
    assert_error!(
      result,
      "step `anc` reads `/pipeline/tmp/pipeline-safety-absent/tt/timetree.nwk`, produced by step `tt`, which is not selected and whose output is absent; include `tt` in --steps or run it first"
    );
  }

  #[test]
  fn test_safety_full_run_has_no_missing_upstream() {
    validate(
      json!({
        "output_all": "tmp/pipeline-safety-absent",
        "steps": [
          { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv" } },
          { "name": "anc", "ancestral": { "tree": "{{ steps.tt.outputs.nwk }}" } }
        ]
      }),
      None,
    )
    .unwrap();
  }

  #[test]
  fn test_safety_step_keep_root_with_reroot_rejected_at_the_key() {
    let config = json!({
      "output_all": "out",
      "steps": [
        { "name": "a", "clock": { "tree": "in.nwk", "metadata": "m.tsv", "keep_root": true, "reroot": "oldest" } }
      ]
    });
    let keep_root = config.to_string().find("\"keep_root\"").unwrap();
    let report = validate(config, None).err().unwrap();
    assert_eq!(
      (
        o!("in pipeline step `a`: invalid configuration: `keep_root` cannot be used together with `reroot`"),
        vec![Some(keep_root)]
      ),
      (
        report_to_string(&report),
        InvalidConfig::problems_of(&report)
          .iter()
          .map(|problem| problem.span.map(|span| span.offset))
          .collect::<Vec<_>>()
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::template_resolves_to_true( true,  Some("in pipeline step `a`: invalid configuration: `keep_root` cannot be used together with `reroot`"))]
  #[case::template_resolves_to_false(false, None)]
  #[trace]
  fn test_safety_step_checks_conflicts_after_interpolation(#[case] keep: bool, #[case] expected: Option<&str>) {
    let result = validate(
      json!({
        "vars": { "keep": keep },
        "output_all": "out",
        "steps": [
          { "name": "a", "clock": { "tree": "in.nwk", "metadata": "m.tsv", "keep_root": "{{ vars.keep }}", "reroot": "oldest" } }
        ]
      }),
      None,
    );
    assert_eq!(
      expected.map(ToOwned::to_owned),
      result.err().map(|report| report_to_string(&report))
    );
  }
}
