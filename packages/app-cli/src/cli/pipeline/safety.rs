use crate::cli::pipeline::inputs::step_inputs;
use crate::cli::pipeline::resolve::{ResolvedPipeline, ResolvedStep};
use app_commands::command_config::CommandConfig;
use app_commands::commands::ancestral::aa_node_data::{
  cds_names, cds_output_paths, cds_selection_names, translation_path,
};
use app_commands::commands::shared::run_paths::RunPath;
use app_commands::config::cli_rules::cli_rule_diagnostics;
use app_commands::config::properties::PathRole;
use app_commands::config::source::{ConfigSource, render_and_bail};
use app_output::output_plan::PathPlaceholder;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use treetime_utils::io::file::is_path_stdout;
use treetime_utils::make_error;

pub(crate) fn validate_plan(pipeline: &ResolvedPipeline, selected: Option<&BTreeSet<String>>) -> Result<(), Report> {
  reject_stdio_paths(pipeline)?;
  let steps: Vec<PlannedStep<'_>> = pipeline
    .steps
    .iter()
    .map(|step| PlannedStep::new(step, selected))
    .try_collect()?;
  reject_output_collisions(&steps)?;
  reject_missing_upstream(&steps, selected)?;
  Ok(())
}

pub(crate) fn reject_step_cli_rule_violations(
  position: usize,
  name: &str,
  command: &CommandConfig,
  source: &ConfigSource,
) -> Result<(), Report> {
  let tag: &str = command.command().into();
  let diags = cli_rule_diagnostics(command.command(), &command.settings()?)?
    .into_iter()
    .map(|mut diag| {
      diag.pointer = diag.pointer.map(|pointer| format!("/steps/{position}/{tag}{pointer}"));
      diag
    })
    .collect();
  render_and_bail(source, "invalid configuration", diags).wrap_err_with(|| format!("in pipeline step `{name}`"))
}

fn reject_stdio_paths(pipeline: &ResolvedPipeline) -> Result<(), Report> {
  for step in &pipeline.steps {
    if step_inputs(&step.command)?.iter().any(RunPath::is_stdio) {
      return make_error!(
        "step `{}` uses `-` (stdin) for an input; a pipeline runs many steps in one process, so steps need real file paths",
        step.name
      );
    }
    if step.outputs.by_selection.values().flatten().any(is_path_stdout) {
      return make_error!(
        "step `{}` uses `-` (stdout) for an output; a pipeline runs many steps in one process, so steps need real file paths",
        step.name
      );
    }
  }
  Ok(())
}

fn reject_output_collisions(steps: &[PlannedStep<'_>]) -> Result<(), Report> {
  let mut owner: BTreeMap<&str, &str> = BTreeMap::new();
  for step in steps {
    for path in &step.outputs {
      if let Some(previous) = owner.insert(path, step.name) {
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

fn reject_missing_upstream(steps: &[PlannedStep<'_>], selected: Option<&BTreeSet<String>>) -> Result<(), Report> {
  let Some(selected) = selected else {
    return Ok(());
  };

  let producers: BTreeMap<&str, &str> = steps
    .iter()
    .flat_map(|step| step.outputs.iter().map(|path| (path.as_str(), step.name)))
    .collect();

  for step in steps.iter().filter(|step| selected.contains(step.name)) {
    for input in &step.inputs {
      if let Some(producer) = producers.get(input.as_str()) {
        if !selected.contains(*producer) && !Path::new(input).exists() {
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

struct PlannedStep<'a> {
  name: &'a str,
  inputs: Vec<String>,
  outputs: Vec<String>,
}

impl<'a> PlannedStep<'a> {
  fn new(step: &'a ResolvedStep, selected: Option<&BTreeSet<String>>) -> Result<Self, Report> {
    let cdses = step_cds_names(step, selected).wrap_err_with(|| format!("in pipeline step `{}`", step.name))?;
    let inputs = step_inputs(&step.command)?
      .into_iter()
      .flat_map(|input| input_paths(input, &cdses))
      .collect();
    let outputs = step
      .outputs
      .by_selection
      .values()
      .flatten()
      .map(|path| output_paths(path, &cdses))
      .flatten_ok()
      .try_collect()
      .wrap_err_with(|| format!("in pipeline step `{}`", step.name))?;
    Ok(Self {
      name: &step.name,
      inputs,
      outputs,
    })
  }
}

fn step_cds_names(step: &ResolvedStep, selected: Option<&BTreeSet<String>>) -> Result<Vec<String>, Report> {
  let settings = step.command.settings()?;
  if selected.is_none_or(|selected| selected.contains(&step.name)) {
    cds_names(&settings)
  } else {
    Ok(cds_selection_names(&settings))
  }
}

fn input_paths(input: RunPath, cdses: &[String]) -> Vec<String> {
  if input.role == PathRole::InputTemplate {
    cdses
      .iter()
      .map(|cds| translation_path(&input.path, cds).to_string_lossy().into_owned())
      .collect()
  } else {
    vec![input.path]
  }
}

fn output_paths(path: &Path, cdses: &[String]) -> Result<Vec<String>, Report> {
  if !PathPlaceholder::Cds.is_in_file_name(path) {
    return Ok(vec![path.to_string_lossy().into_owned()]);
  }
  Ok(
    cds_output_paths(path, cdses)?
      .values()
      .map(|path| path.to_string_lossy().into_owned())
      .collect(),
  )
}

#[cfg(test)]
mod tests {
  const PIPELINE_FOLDER: &str = "/pipeline";

  use crate::cli::pipeline::safety::tests::helpers::{ancestral_step, validate, validate_in, write_annotation};
  use app_commands::config::source::InvalidConfig;
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use tempfile::tempdir;
  use treetime_utils::error::report_to_string;
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_safety_output_collision_rejected() {
    let result = validate(
      json!({
        "steps": [
          { "name": "a", "ancestral": { "alignment": ["a.fasta"], "tree": "in.nwk", "output_tree_nwk": "shared.nwk" } },
          { "name": "b", "ancestral": { "alignment": ["a.fasta"], "tree": "in.nwk", "output_tree_nwk": "shared.nwk" } }
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
  fn test_safety_output_collision_of_a_per_alignment_file_rejected() {
    let result = validate(
      json!({
        "steps": [
          { "name": "a", "ancestral": { "tree": "in.nwk", "alignment": ["x.fasta", "y.fasta"], "output_sequences": "seq.{alignment}.fasta" } },
          { "name": "b", "ancestral": { "alignment": ["a.fasta"], "tree": "in.nwk", "output_sequences": "seq.2.fasta" } }
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
          { "name": "a", "ancestral": { "alignment": ["a.fasta"], "tree": "same.nwk", "output_tree_nwk": "same.nwk" } }
        ]
      }),
      None,
    );
    assert_error!(
      result,
      "in pipeline step `a`: --tree '/pipeline/same.nwk' names the same file as the output --output-tree-nwk \
       '/pipeline/same.nwk': the run would overwrite its own input"
    );
  }

  #[test]
  fn test_safety_translation_template_that_matches_an_output_rejected() {
    let result = validate(
      json!({
        "steps": [
          { "name": "a", "ancestral": {
            "tree": "in.nwk",
            "alignment": ["aln.fasta"],
            "translations": "aa/{cds}.fasta",
            "cds_selection": ["S"],
            "output_sequences": "aa/S.fasta"
          } }
        ]
      }),
      None,
    );
    assert_error!(
      result,
      "in pipeline step `a`: --translations '/pipeline/aa/S.fasta' names the same file as the output \
       --output-sequences '/pipeline/aa/S.fasta': the run would overwrite its own input"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::equal_translation_templates(
    json!({ "translations": "tr/{cds}.fasta", "cds_selection": ["F"], "output_translations": "aa/{cds}.fasta" }),
    json!({ "translations": "tr/{cds}.fasta", "cds_selection": ["F"], "output_translations": "aa/{cds}.fasta" }),
    "aa/F.fasta",
  )]
  #[case::translation_template_expands_to_a_plain_output(
    json!({ "translations": "tr/{cds}.fasta", "cds_selection": ["F"], "output_translations": "aa/{cds}.fasta" }),
    json!({ "output_sequences": "aa/F.fasta" }),
    "aa/F.fasta",
  )]
  #[case::cds_names_of_the_annotation(
    json!({ "translations": "tr/{cds}.fasta", "annotation": "genes.gff", "output_translations": "aa/{cds}.fasta" }),
    json!({ "output_sequences": "aa/G.fasta" }),
    "aa/G.fasta",
  )]
  #[trace]
  fn test_safety_translation_output_collision_rejected(#[case] a: Value, #[case] b: Value, #[case] shared: &str) {
    let dir = tempdir().unwrap();
    write_annotation(dir.path(), "genes.gff", &["F", "G"]);
    let result = validate_in(
      dir.path(),
      json!({ "steps": [ancestral_step("a", a), ancestral_step("b", b)] }),
      None,
    );
    assert_error!(
      result,
      format!(
        "steps `a` and `b` both write `{}`; each output path must be written by only one step",
        dir.path().join(shared).display()
      )
    );
  }

  #[test]
  fn test_safety_missing_upstream_translation_rejected() {
    let dir = tempdir().unwrap();
    let result = validate_in(
      dir.path(),
      json!({
        "steps": [
          ancestral_step("a", json!({ "translations": "tr/{cds}.fasta", "cds_selection": ["F"], "output_translations": "aa/{cds}.fasta" })),
          ancestral_step("b", json!({ "translations": "aa/{cds}.fasta", "cds_selection": ["F"], "output_sequences": "b.fasta" }))
        ]
      }),
      Some(&btreeset! { o!("b") }),
    );
    assert_error!(
      result,
      format!(
        "step `b` reads `{}`, produced by step `a`, which is not selected and whose output is absent; \
         include `a` in --steps or run it first",
        dir.path().join("aa/F.fasta").display()
      )
    );
  }

  #[test]
  fn test_safety_annotation_of_an_unselected_step_is_not_read() {
    let dir = tempdir().unwrap();
    validate_in(
      dir.path(),
      json!({
        "steps": [
          ancestral_step("a", json!({ "output_sequences": "a.fasta" })),
          ancestral_step("b", json!({ "translations": "tr/{cds}.fasta", "annotation": "absent.gff", "output_translations": "aa/{cds}.fasta" }))
        ]
      }),
      Some(&btreeset! { o!("a") }),
    )
    .unwrap();
  }

  #[test]
  fn test_safety_absent_annotation_of_a_selected_step_rejected() {
    let dir = tempdir().unwrap();
    let result = validate_in(
      dir.path(),
      json!({
        "steps": [
          ancestral_step("b", json!({ "translations": "tr/{cds}.fasta", "annotation": "absent.gff", "output_translations": "aa/{cds}.fasta" }))
        ]
      }),
      None,
    );
    assert_error!(
      result,
      format!(
        "in pipeline step `b`: When opening file '{}': No such file or directory (os error 2)",
        dir.path().join("absent.gff").display()
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tree(        json!({ "ancestral": { "alignment": ["a.fasta"], "tree": "-", "output_model": "g.json" } }))]
  #[case::remove_file( json!({ "prune": { "tree": "in.nwk", "remove_file": "-", "output_tree_nwk": "out.nwk" } }))]
  #[trace]
  fn test_safety_stdio_input_rejected(#[case] step: Value) {
    let mut step = step;
    step["name"] = json!("a");
    let result = validate(json!({ "steps": [step] }), None);
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
          { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv", "sequence_length": 1000 } },
          { "name": "anc", "ancestral": { "alignment": ["a.fasta"], "tree": "{{ steps.tt.outputs.nwk }}" } }
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
          { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv", "sequence_length": 1000 } },
          { "name": "anc", "ancestral": { "alignment": ["a.fasta"], "tree": "{{ steps.tt.outputs.nwk }}" } }
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

  #[test]
  fn test_safety_step_partition_with_model_params_reports_the_conflict() {
    let result = validate(
      json!({
        "output_all": "out",
        "steps": [
          { "name": "a", "ancestral": { "tree": "in.nwk", "alignment": ["a.fasta"], "partition": [{}], "model_params": ["mu=2"] } }
        ]
      }),
      None,
    );
    assert_error!(
      result,
      "in pipeline step `a`: invalid configuration: `model_params` cannot be used together with `partition`"
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

  mod helpers {
    use crate::cli::pipeline::resolve::{PipelineDoc, resolve_pipeline};
    use crate::cli::pipeline::safety::validate_plan;
    use app_commands::config::source::ConfigSource;
    use eyre::Report;
    use itertools::Itertools;
    use serde_json::{Value, json};
    use std::collections::BTreeSet;
    use std::fs;
    use std::path::Path;

    pub(super) fn validate(config: Value, selected: Option<&BTreeSet<String>>) -> Result<(), Report> {
      validate_in(Path::new(super::PIPELINE_FOLDER), config, selected)
    }

    pub(super) fn validate_in(folder: &Path, config: Value, selected: Option<&BTreeSet<String>>) -> Result<(), Report> {
      let source = ConfigSource::new("pipeline.json", config.to_string());
      let doc = PipelineDoc::from_value(config).unwrap();
      let pipeline = resolve_pipeline(&doc, &source, &json!({}), folder, None)?;
      validate_plan(&pipeline, selected)
    }

    pub(super) fn ancestral_step(name: &str, settings: Value) -> Value {
      let Value::Object(settings) = settings else {
        panic!("the settings of a step must be an object, got {settings}")
      };
      let mut ancestral = json!({ "alignment": ["aln.fasta"], "tree": "in.nwk" });
      ancestral.as_object_mut().unwrap().extend(settings);
      json!({ "name": name, "ancestral": ancestral })
    }

    pub(super) fn write_annotation(folder: &Path, name: &str, cdses: &[&str]) {
      let features = cdses
        .iter()
        .map(|cds| format!("ref\t.\tCDS\t1\t9\t.\t+\t0\tName={cds}"))
        .join("\n");
      fs::write(folder.join(name), format!("##gff-version 3\n{features}\n")).unwrap();
    }
  }
}
