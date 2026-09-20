use crate::cli::pipeline::interpolate::{Interpolator, map_string_leaves, resolve_vars, template_context};
use crate::cli::pipeline::suggest::{suggestion_suffix, valid_values};
use crate::cli::pipeline::types::{PipelineStepCommand, RawStep};
use app_output::output_plan::OutputSelection;
use eyre::Report;
use itertools::Itertools;
use regex::{Regex, regex};
use serde_json::{Map, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use treetime_utils::{make_error, make_report};

pub(crate) const TOP_LEVEL_KEYS: [&str; 4] = ["$schema", "vars", "output_all", "steps"];

pub(crate) const CDS_PLACEHOLDERS: [&str; 2] = ["{cds}", "%GENE"];

pub(crate) fn step_ref() -> &'static Regex {
  regex!(
    r#"\{\{\s*steps\.([A-Za-z0-9_]+)\.(?:outputs\.([A-Za-z0-9_\-]+)|outputs\[\s*['"]([A-Za-z0-9_\-]+)['"]\s*\]|(output_all))\s*\}\}"#
  )
}

pub struct PipelineDoc {
  pub schema_ref: Option<String>,
  pub vars: Map<String, Value>,
  pub output_all: Option<String>,
  pub steps: Vec<RawStep>,
}

impl PipelineDoc {
  pub fn from_value(value: Value) -> Result<Self, Report> {
    let Value::Object(mut map) = value else {
      return make_error!("a pipeline config must be a mapping with `steps`");
    };

    for key in map.keys() {
      if !TOP_LEVEL_KEYS.contains(&key.as_str()) {
        return make_error!(
          "unknown top-level key `{key}`; {}",
          suggestion_suffix(key, &TOP_LEVEL_KEYS)
        );
      }
    }

    let schema_ref = match map.remove("$schema") {
      Some(Value::String(reference)) => Some(reference),
      Some(_) => return make_error!("top-level `$schema` must be a string"),
      None => None,
    };

    let vars = match map.remove("vars") {
      Some(Value::Object(vars)) => vars,
      Some(_) => return make_error!("top-level `vars` must be a mapping"),
      None => Map::new(),
    };

    let output_all = match map.remove("output_all") {
      Some(Value::String(dir)) => Some(dir),
      Some(_) => return make_error!("top-level `output_all` must be a string"),
      None => None,
    };

    let steps_value = map
      .remove("steps")
      .ok_or_else(|| make_report!("pipeline config has no `steps`"))?;
    let Value::Array(step_values) = steps_value else {
      return make_error!("`steps` must be a list of steps");
    };

    let steps: Vec<RawStep> = step_values.into_iter().map(RawStep::from_value).try_collect()?;

    let mut seen = BTreeSet::new();
    for step in &steps {
      if !seen.insert(step.name.clone()) {
        return make_error!("duplicate step name `{}`; step names must be unique", step.name);
      }
    }

    Ok(Self {
      schema_ref,
      vars,
      output_all,
      steps,
    })
  }
}

pub struct StepOutputs {
  pub output_all: Option<PathBuf>,
  pub by_selection: BTreeMap<OutputSelection, Vec<PathBuf>>,
}

pub struct ResolvedStep {
  pub name: String,
  pub command: PipelineStepCommand,
  pub outputs: StepOutputs,
}

pub struct ResolvedPipeline {
  pub workdir: Option<PathBuf>,
  pub steps: Vec<ResolvedStep>,
}

pub fn resolve_pipeline(doc: &PipelineDoc, env: &Value) -> Result<ResolvedPipeline, Report> {
  let interp = Interpolator::default();
  let vars = resolve_vars(&interp, &doc.vars, env)?;

  let value_context = template_context(&vars, env);
  let workdir = match &doc.output_all {
    Some(template) => Some(interpolate_to_path(&interp, template, &value_context, "output_all")?),
    None => None,
  };

  let all_names: BTreeSet<&str> = doc.steps.iter().map(|step| step.name.as_str()).collect();
  let mut resolved_steps: Vec<ResolvedStep> = Vec::new();
  let mut index: BTreeMap<String, usize> = BTreeMap::new();

  for raw in &doc.steps {
    let payload = substitute_step_refs(&raw.payload, &resolved_steps, &index, &raw.name, &all_names)?;
    let mut payload = interp
      .interpolate_value(&payload, &value_context)
      .map_err(|err| step_error(&raw.name, &err))?;

    if let Some(dir) = &workdir {
      set_output_all_if_absent(&mut payload, &dir.join(&raw.name));
    }

    let command =
      PipelineStepCommand::from_tag_and_value(&raw.tag, payload).map_err(|err| step_error(&raw.name, &err))?;
    let resolved = command.resolve_outputs().map_err(|err| step_error(&raw.name, &err))?;
    let outputs = StepOutputs {
      output_all: command.output_all().map(Path::to_path_buf),
      by_selection: resolved.paths_by_selection(),
    };

    index.insert(raw.name.clone(), resolved_steps.len());
    resolved_steps.push(ResolvedStep {
      name: raw.name.clone(),
      command,
      outputs,
    });
  }

  Ok(ResolvedPipeline {
    workdir,
    steps: resolved_steps,
  })
}

fn interpolate_to_path(interp: &Interpolator, template: &str, context: &Value, field: &str) -> Result<PathBuf, Report> {
  match interp.interpolate_str(template, context)? {
    Value::String(path) => Ok(PathBuf::from(path)),
    other => make_error!("`{field}` must resolve to a string path, got {other}"),
  }
}

fn substitute_step_refs(
  value: &Value,
  resolved: &[ResolvedStep],
  index: &BTreeMap<String, usize>,
  step: &str,
  all_names: &BTreeSet<&str>,
) -> Result<Value, Report> {
  map_string_leaves(value, &mut |leaf| {
    Ok(Value::String(substitute_step_refs_in_leaf(
      leaf, resolved, index, step, all_names,
    )?))
  })
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
fn substitute_step_refs_in_leaf(
  leaf: &str,
  resolved: &[ResolvedStep],
  index: &BTreeMap<String, usize>,
  step: &str,
  all_names: &BTreeSet<&str>,
) -> Result<String, Report> {
  let mut result = String::new();
  let mut last = 0;
  for captures in step_ref().captures_iter(leaf) {
    let whole = captures.get(0).expect("regex match has group 0");
    result.push_str(leaf.get(last..whole.start()).unwrap_or_default());

    let producer = resolve_producer(&captures[1], resolved, index, step, all_names)?;
    let path = if captures.get(4).is_some() {
      producer_output_all(step, producer)?
    } else {
      let selection = captures
        .get(2)
        .or_else(|| captures.get(3))
        .expect("regex has a selection group")
        .as_str();
      resolve_selection_path(step, producer, selection)?
    };
    result.push_str(&path);
    last = whole.end();
  }
  result.push_str(leaf.get(last..).unwrap_or_default());
  Ok(result)
}

fn resolve_producer<'a>(
  producer: &str,
  resolved: &'a [ResolvedStep],
  index: &BTreeMap<String, usize>,
  step: &str,
  all_names: &BTreeSet<&str>,
) -> Result<&'a ResolvedStep, Report> {
  if let Some(&position) = index.get(producer) {
    return Ok(&resolved[position]);
  }
  if producer == step || all_names.contains(producer) {
    return make_error!(
      "step `{step}` references `{producer}`, which is not an earlier step; steps may only reference steps before them"
    );
  }
  let earlier: Vec<&str> = index.keys().map(String::as_str).collect();
  make_error!(
    "step `{step}` references unknown step `{producer}`; {}",
    suggestion_suffix(producer, &earlier)
  )
}

fn producer_output_all(step: &str, producer: &ResolvedStep) -> Result<String, Report> {
  match &producer.outputs.output_all {
    Some(dir) => Ok(dir.to_string_lossy().into_owned()),
    None => make_error!(
      "step `{step}` references `steps.{}.output_all`, but that step has no output directory",
      producer.name
    ),
  }
}

fn resolve_selection_path(step: &str, producer: &ResolvedStep, selection: &str) -> Result<String, Report> {
  let by_selection = &producer.outputs.by_selection;
  let produced_tags: Vec<String> = by_selection.keys().map(|sel| selection_tag(*sel)).collect();
  let produced_tags: Vec<&str> = produced_tags.iter().map(String::as_str).collect();

  let Some(parsed) = parse_selection(selection) else {
    return make_error!(
      "step `{step}` references `steps.{}.outputs.{selection}`, an unknown output; {}",
      producer.name,
      suggestion_suffix(selection, &produced_tags)
    );
  };

  match by_selection.get(&parsed).map(Vec::as_slice) {
    None | Some([]) => make_error!(
      "step `{step}` references `{selection}` from step `{}`, which does not produce it; it produces: {}",
      producer.name,
      valid_values(&produced_tags)
    ),
    Some([path]) if is_template_path(path) => make_error!(
      "step `{step}` references `{selection}` from step `{}`, which resolves to a per-CDS template (`{}`); \
       reference a specific CDS output path instead",
      producer.name,
      path.display()
    ),
    Some([path]) => Ok(path.to_string_lossy().into_owned()),
    Some(paths) => make_error!(
      "step `{step}` references `{selection}` from step `{}`, which resolves to more than one file ({}); \
       pin a single style or CDS, or reference a specific path",
      producer.name,
      paths.iter().map(|path| format!("`{}`", path.display())).join(", ")
    ),
  }
}

fn set_output_all_if_absent(payload: &mut Value, dir: &Path) {
  let Value::Object(map) = payload else {
    return;
  };
  let already_set = matches!(map.get("output_all"), Some(value) if !value.is_null());
  if !already_set {
    map.insert(
      "output_all".to_owned(),
      Value::String(dir.to_string_lossy().into_owned()),
    );
  }
}

fn selection_tag(selection: OutputSelection) -> String {
  serde_json::to_value(selection)
    .ok()
    .and_then(|value| value.as_str().map(str::to_owned))
    .unwrap_or_default()
}

fn parse_selection(tag: &str) -> Option<OutputSelection> {
  serde_json::from_value(Value::String(tag.to_owned())).ok()
}

pub(crate) fn is_template_path(path: &Path) -> bool {
  let path = path.to_string_lossy();
  CDS_PLACEHOLDERS.iter().any(|placeholder| path.contains(placeholder))
}

fn step_error(step: &str, err: &Report) -> Report {
  make_report!("in pipeline step `{step}`: {err}")
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use treetime_utils::assert_error;

  fn resolve(config: Value) -> Result<ResolvedPipeline, Report> {
    let doc = PipelineDoc::from_value(config)?;
    resolve_pipeline(&doc, &json!({}))
  }

  fn tree_of(step: &ResolvedStep) -> String {
    match &step.command {
      PipelineStepCommand::Ancestral(args) => args
        .tree
        .as_deref()
        .expect("ancestral step has a tree")
        .to_string_lossy()
        .into_owned(),
      _ => panic!("expected ancestral step"),
    }
  }

  fn step<'a>(pipeline: &'a ResolvedPipeline, name: &str) -> &'a ResolvedStep {
    pipeline
      .steps
      .iter()
      .find(|step| step.name == name)
      .expect("step present")
  }

  #[test]
  fn test_resolve_chain_single_file_reference() {
    let config = json!({
      "output_all": "tmp/run",
      "steps": [
        { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv" } },
        { "name": "anc", "ancestral": { "tree": "{{ steps.tt.outputs.nwk }}", "method_anc": "marginal" } }
      ]
    });
    let resolved = resolve(config).unwrap();
    assert_eq!("tmp/run/tt/timetree.nwk", tree_of(step(&resolved, "anc")));
  }

  #[test]
  fn test_resolve_chain_uses_explicit_step_output_all() {
    let config = json!({
      "steps": [
        { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv", "output_all": "out/tt" } },
        { "name": "anc", "ancestral": { "tree": "{{ steps.tt.outputs.nwk }}", "output_all": "out/anc" } }
      ]
    });
    let resolved = resolve(config).unwrap();
    assert_eq!("out/tt/timetree.nwk", tree_of(step(&resolved, "anc")));
  }

  #[test]
  fn test_resolve_chain_multi_style_reference_is_ambiguous() {
    let config = json!({
      "output_all": "tmp/run",
      "steps": [
        { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv", "output_nwk_style": ["plain", "beast"] } },
        { "name": "anc", "ancestral": { "tree": "{{ steps.tt.outputs.nwk }}" } }
      ]
    });
    let result = resolve(config);
    assert_error!(
      result,
      "step `anc` references `nwk` from step `tt`, which resolves to more than one file (`tmp/run/tt/timetree.annotated.nwk`, `tmp/run/tt/timetree.nwk`); pin a single style or CDS, or reference a specific path"
    );
  }

  #[test]
  fn test_resolve_chain_per_cds_template_reference_rejected() {
    let config = json!({
      "output_all": "tmp/run",
      "steps": [
        { "name": "anc", "ancestral": { "tree": "in.nwk", "output_reconstructed_aa_fasta": "out/{cds}.fasta" } },
        { "name": "anc2", "ancestral": { "tree": "{{ steps.anc.outputs.reconstructed-aa-fasta }}" } }
      ]
    });
    let result = resolve(config);
    assert_error!(
      result,
      "step `anc2` references `reconstructed-aa-fasta` from step `anc`, which resolves to a per-CDS template (`out/{cds}.fasta`); reference a specific CDS output path instead"
    );
  }

  #[test]
  fn test_resolve_chain_forward_reference_rejected() {
    let config = json!({
      "output_all": "tmp/run",
      "steps": [
        { "name": "anc", "ancestral": { "tree": "{{ steps.tt.outputs.nwk }}" } },
        { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv" } }
      ]
    });
    let result = resolve(config);
    assert_error!(
      result,
      "step `anc` references `tt`, which is not an earlier step; steps may only reference steps before them"
    );
  }

  #[test]
  fn test_resolve_chain_unknown_step_reference_suggests() {
    let config = json!({
      "output_all": "tmp/run",
      "steps": [
        { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv" } },
        { "name": "anc", "ancestral": { "tree": "{{ steps.tta.outputs.nwk }}" } }
      ]
    });
    let result = resolve(config);
    assert_error!(
      result,
      "step `anc` references unknown step `tta`; did you mean `tt`? Valid values: `tt`"
    );
  }

  #[test]
  fn test_resolve_chain_unproduced_selection_rejected() {
    let config = json!({
      "output_all": "tmp/run",
      "steps": [
        { "name": "tt", "timetree": { "tree": "in.nwk", "metadata": "m.tsv" } },
        { "name": "anc", "ancestral": { "tree": "{{ steps.tt.outputs.traits-csv }}" } }
      ]
    });
    let result = resolve(config);
    assert_error!(
      result,
      "step `anc` references `traits-csv` from step `tt`, which does not produce it; it produces: `augur-node-data`, `auspice`, `clock-model`, `coalescent-tsv`, `gtr`, `nexus`, `nwk`, `reconstructed-nuc-fasta`"
    );
  }

  #[test]
  fn test_resolve_duplicate_step_names_rejected() {
    let config = json!({
      "steps": [
        { "name": "tt", "timetree": { "tree": "a.nwk" } },
        { "name": "tt", "timetree": { "tree": "b.nwk" } }
      ]
    });
    let result = PipelineDoc::from_value(config);
    assert_error!(result, "duplicate step name `tt`; step names must be unique");
  }

  #[test]
  fn test_resolve_unknown_top_level_key_rejected() {
    let config = json!({ "step": [] });
    let result = PipelineDoc::from_value(config);
    assert_error!(
      result,
      "unknown top-level key `step`; did you mean `steps`? Valid values: `$schema`, `output_all`, `steps`, `vars`"
    );
  }
}
