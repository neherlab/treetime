use crate::cli::config::config_folder;
use crate::cli::diagnostics::entry::check_pipeline;
use crate::cli::pipeline::check::print_pipeline_plan;
use crate::cli::pipeline::resolve::{PipelineDoc, ResolvedPipeline, ResolvedStep, resolve_pipeline};
use crate::cli::pipeline::safety::validate_plan;
use crate::cli::treetime_cli::TreetimePipelineArgs;
use app_commands::config::source::{ConfigSource, parse_config_document};
use app_commands::config::suggest::suggestion_suffix;
use eyre::Report;
use itertools::Itertools;
use serde_json::{Map, Value};
use std::collections::BTreeSet;
use std::path::Path;
use treetime::cancel::NoopCancel;
use treetime::progress::{LogSink, StageSink};
use treetime_utils::io::fs::{absolute_path, read_file_to_string};
use treetime_utils::make_error;

pub(crate) fn run_pipeline_command(
  args: &TreetimePipelineArgs,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let output_all = args.output_all.as_deref().map(absolute_path).transpose()?;
  let pipeline = load_pipeline(&args.config, output_all.as_deref())?;
  let selected = (!args.steps.is_empty()).then(|| args.steps.iter().cloned().collect::<BTreeSet<String>>());
  validate_plan(&pipeline, selected.as_ref())?;
  if args.check {
    print_pipeline_plan(&pipeline, selected.as_ref())
  } else {
    run_pipeline(&pipeline, selected.as_ref(), stages, log)
  }
}

pub(crate) fn load_pipeline(config: &Path, output_all: Option<&Path>) -> Result<ResolvedPipeline, Report> {
  let text = read_file_to_string(config)?;
  let source = ConfigSource::new(config.display().to_string(), text.clone());

  let value = parse_config_document(&source, &text)?;

  check_pipeline(&source, &value)?;
  let doc = PipelineDoc::from_value(value)?;
  resolve_pipeline(&doc, &process_env(), &config_folder(config)?, output_all)
}

fn process_env() -> Value {
  let env = std::env::vars()
    .map(|(key, value)| (key, Value::String(value)))
    .collect::<Map<_, _>>();
  Value::Object(env)
}

pub(crate) fn run_pipeline(
  pipeline: &ResolvedPipeline,
  selected: Option<&BTreeSet<String>>,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let steps = select_steps(pipeline, selected)?;

  let mut completed: Vec<&ResolvedStep> = Vec::new();
  for (position, step) in steps.iter().enumerate() {
    if let Err(err) = run_step(step, stages, log) {
      let remaining = steps[position..].iter().map(|step| step.name.as_str()).join(",");
      return Err(err.wrap_err(failure_report(step, &completed, &remaining)));
    }
    completed.push(step);
  }

  Ok(())
}

pub(crate) fn select_steps<'a>(
  pipeline: &'a ResolvedPipeline,
  selected: Option<&BTreeSet<String>>,
) -> Result<Vec<&'a ResolvedStep>, Report> {
  let Some(selected) = selected else {
    return Ok(pipeline.steps.iter().collect());
  };

  let names: BTreeSet<&str> = pipeline.steps.iter().map(|step| step.name.as_str()).collect();
  for name in selected {
    if !names.contains(name.as_str()) {
      let candidates: Vec<&str> = names.iter().copied().collect();
      return make_error!(
        "unknown step `{name}` in --steps; {}",
        suggestion_suffix(name, &candidates)
      );
    }
  }

  Ok(
    pipeline
      .steps
      .iter()
      .filter(|step| selected.contains(&step.name))
      .collect(),
  )
}

fn run_step(step: &ResolvedStep, stages: &dyn StageSink, log: &dyn LogSink) -> Result<(), Report> {
  step.command.args()?.execute(&NoopCancel, stages, log)
}

fn failure_report(failed: &ResolvedStep, completed: &[&ResolvedStep], remaining: &str) -> String {
  let done = if completed.is_empty() {
    "none".to_owned()
  } else {
    completed
      .iter()
      .map(|step| format!("`{}`{}", step.name, output_dir_hint(step.outputs.output_all.as_deref())))
      .join(", ")
  };
  format!(
    "pipeline step `{}` failed; completed steps: {done}; resume the remaining steps with --steps={remaining}",
    failed.name
  )
}

fn output_dir_hint(output_all: Option<&Path>) -> String {
  output_all.map_or_else(String::new, |dir| format!(" (outputs in {})", dir.display()))
}

#[cfg(test)]
mod tests {
  const PIPELINE_FOLDER: &str = "/pipeline";

  use super::*;
  use crate::cli::pipeline::resolve::PipelineDoc;
  use indoc::indoc;
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::path::PathBuf;
  use std::{env, fs, iter};
  use tempfile::tempdir;
  use treetime_utils::assert_error;

  fn resolve(config: Value) -> ResolvedPipeline {
    let doc = PipelineDoc::from_value(config).unwrap();
    resolve_pipeline(&doc, &json!({}), Path::new(PIPELINE_FOLDER), None).unwrap()
  }

  fn config() -> Value {
    json!({
      "output_all": "tmp/run",
      "steps": [
        { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv" } },
        { "name": "anc", "ancestral": { "tree": "{{ steps.tt.outputs.nwk }}" } }
      ]
    })
  }

  #[test]
  fn test_runner_select_all_steps_in_order() {
    let pipeline = resolve(config());
    let selected = select_steps(&pipeline, None).unwrap();
    assert_eq!(
      vec!["tt", "anc"],
      selected.iter().map(|step| step.name.as_str()).collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_runner_select_subset_keeps_list_order() {
    let pipeline = resolve(config());
    let selected = select_steps(&pipeline, Some(&btreeset! { "anc".to_owned() })).unwrap();
    assert_eq!(
      vec!["anc"],
      selected.iter().map(|step| step.name.as_str()).collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_runner_select_unknown_step_errors() {
    let pipeline = resolve(config());
    let result = select_steps(&pipeline, Some(&btreeset! { "tta".to_owned() }));
    assert_error!(
      result,
      "unknown step `tta` in --steps; did you mean `tt`? Valid values: `anc`, `tt`"
    );
  }

  #[test]
  fn test_runner_load_pipeline_resolves_a_relative_config_and_a_relative_output_folder() {
    let dir = tempdir().unwrap();
    let config = dir.path().join("pipeline.yaml");
    fs::write(
      &config,
      indoc! {r"
        steps:
          - name: tt
            timetree: { tree: tree.nwk, metadata: metadata.tsv }
          - name: anc
            ancestral: { tree: '{{ steps.tt.outputs.nwk }}', alignment: [aln.fasta.xz] }
      "},
    )
    .unwrap();
    let cwd = env::current_dir().unwrap();
    let relative_config = relative_from(&cwd, &config);
    let relative_out = relative_from(&cwd, &dir.path().join("out"));
    let output_all = absolute_path(&relative_out).unwrap();

    let pipeline = load_pipeline(&relative_config, Some(&output_all)).unwrap();

    let anc = pipeline.steps.iter().find(|step| step.name == "anc").unwrap();
    let args = Value::Object(anc.command.settings().unwrap());
    let config_folder = cwd.join(&relative_config).parent().unwrap().to_path_buf();
    assert_eq!(
      (
        json!(cwd.join(&relative_out).join("tt/timetree.nwk")),
        json!([config_folder.join("aln.fasta.xz")]),
        json!(cwd.join(&relative_out).join("anc"))
      ),
      (
        args["tree"].clone(),
        args["alignment"].clone(),
        args["output_all"].clone()
      )
    );
  }

  fn relative_from(cwd: &Path, target: &Path) -> PathBuf {
    let up = cwd.components().count() - 1;
    let target = target.strip_prefix("/").unwrap();
    iter::repeat_n(Path::new(".."), up).collect::<PathBuf>().join(target)
  }
}
