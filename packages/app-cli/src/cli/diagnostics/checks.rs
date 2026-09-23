use crate::cli::diagnostics::source::{RawDiagnostic, escape_pointer};
use crate::cli::pipeline::interpolate::{Interpolator, NAMESPACES};
use crate::cli::pipeline::resolve::{TOP_LEVEL_KEYS, step_ref};
use crate::cli::pipeline::suggest::suggestion_suffix;
use crate::cli::pipeline::types::{COMMAND_TAGS, SCHEMA_KEY, commands_list};
use crate::cli::schema::command_schema_for;
use itertools::Itertools;
use jsonschema::error::ValidationErrorKind;
use schemars::Schema;
use serde_json::{Map, Value};
use std::collections::{BTreeMap, BTreeSet};

pub fn pipeline_schema_diagnostics(value: &Value) -> Vec<RawDiagnostic> {
  let mut diags = Vec::new();
  let Some(steps) = value.get("steps").and_then(Value::as_array) else {
    return diags;
  };
  for (position, step) in steps.iter().enumerate() {
    let Some(map) = step.as_object() else {
      continue;
    };
    let tags: Vec<&String> = map
      .keys()
      .filter(|key| key.as_str() != "name" && key.as_str() != SCHEMA_KEY)
      .collect();
    let [tag] = tags.as_slice() else {
      continue;
    };
    let Some(schema) = command_schema_for(tag) else {
      continue;
    };
    let base = format!("/steps/{position}/{tag}");
    diags.extend(schema_diagnostics_prefixed(&map[*tag], &schema, true, &base));
  }
  diags
}

pub fn schema_diagnostics(value: &Value, schema: &Schema, skip_templates: bool) -> Vec<RawDiagnostic> {
  schema_diagnostics_prefixed(value, schema, skip_templates, "")
}

fn schema_diagnostics_prefixed(value: &Value, schema: &Schema, skip_templates: bool, base: &str) -> Vec<RawDiagnostic> {
  let schema_value = match serde_json::to_value(schema) {
    Ok(schema_value) => schema_value,
    Err(err) => {
      return vec![RawDiagnostic::new(
        "config::internal",
        format!("could not build schema: {err}"),
      )];
    },
  };
  let validator = match jsonschema::validator_for(&schema_value) {
    Ok(validator) => validator,
    Err(err) => return vec![RawDiagnostic::new("config::internal", format!("invalid schema: {err}"))],
  };

  let mut diags = Vec::new();
  for error in validator.iter_errors(value) {
    if skip_templates {
      if let Value::String(text) = error.instance().as_ref() {
        if text.contains("{{") {
          continue;
        }
      }
    }

    let pointer = format!("{base}{}", error.instance_path().as_str());
    match error.kind() {
      ValidationErrorKind::Enum { options } => {
        let candidates = string_options(options);
        let candidate_refs: Vec<&str> = candidates.iter().map(String::as_str).collect();
        let bad = instance_string(error.instance().as_ref());
        diags.push(
          RawDiagnostic::new("config::enum", format!("`{bad}` is not a valid value"))
            .at(pointer)
            .help(suggestion_suffix(&bad, &candidate_refs)),
        );
      },
      ValidationErrorKind::Type { .. } => {
        diags.push(RawDiagnostic::new("config::type", error.to_string()).at(pointer));
      },
      ValidationErrorKind::Required { property } => {
        let property = property.as_str().unwrap_or_default();
        diags.push(RawDiagnostic::new("config::required", format!("missing required field `{property}`")).at(pointer));
      },
      ValidationErrorKind::AdditionalProperties { unexpected } => {
        for key in unexpected {
          diags.push(
            RawDiagnostic::new("config::unknown-field", format!("unknown field `{key}`"))
              .at(format!("{pointer}/{key}"))
              .key_span(),
          );
        }
      },
      _ => {
        diags.push(RawDiagnostic::new("config::schema", error.to_string()).at(pointer));
      },
    }
  }
  diags
}

pub fn pipeline_structural_diagnostics(value: &Value) -> Vec<RawDiagnostic> {
  let mut diags = Vec::new();
  let Value::Object(map) = value else {
    diags.push(RawDiagnostic::new("config::shape", "a pipeline config must be a mapping with `steps`").at(""));
    return diags;
  };

  for key in map.keys() {
    if !TOP_LEVEL_KEYS.contains(&key.as_str()) {
      diags.push(
        RawDiagnostic::new("config::unknown-field", format!("unknown top-level key `{key}`"))
          .at(format!("/{key}"))
          .key_span()
          .help(suggestion_suffix(key, &TOP_LEVEL_KEYS)),
      );
    }
  }

  if let Some(vars) = map.get("vars") {
    if !vars.is_object() {
      diags.push(RawDiagnostic::new("config::type", "top-level `vars` must be a mapping").at("/vars"));
    }
  }
  if let Some(output_all) = map.get("output_all") {
    if !output_all.is_string() {
      diags.push(RawDiagnostic::new("config::type", "top-level `output_all` must be a string").at("/output_all"));
    }
  }

  match map.get("steps") {
    None => diags.push(RawDiagnostic::new("config::missing-steps", "pipeline config has no `steps`").at("")),
    Some(Value::Array(steps)) => {
      let mut seen = BTreeSet::new();
      for (position, step) in steps.iter().enumerate() {
        structural_step_diagnostics(&mut diags, position, step, &mut seen);
      }
    },
    Some(_) => diags.push(RawDiagnostic::new("config::type", "`steps` must be a list of steps").at("/steps")),
  }

  diags
}

fn structural_step_diagnostics(
  diags: &mut Vec<RawDiagnostic>,
  position: usize,
  step: &Value,
  seen: &mut BTreeSet<String>,
) {
  let base = format!("/steps/{position}");
  let Value::Object(map) = step else {
    diags.push(
      RawDiagnostic::new(
        "config::shape",
        "a pipeline step must be a mapping with a `name` and one command",
      )
      .at(base),
    );
    return;
  };

  let name = match map.get("name") {
    Some(Value::String(name)) => Some(name.clone()),
    Some(_) => {
      diags
        .push(RawDiagnostic::new("config::type", "pipeline step `name` must be a string").at(format!("{base}/name")));
      None
    },
    None => {
      diags.push(RawDiagnostic::new("config::missing-name", "pipeline step is missing a `name`").at(base.clone()));
      None
    },
  };
  if let Some(name) = &name {
    if !seen.insert(name.clone()) {
      diags.push(
        RawDiagnostic::new(
          "config::duplicate-step",
          format!("duplicate step name `{name}`; step names must be unique"),
        )
        .at(format!("{base}/name")),
      );
    }
  }

  let tags: Vec<&String> = map
    .keys()
    .filter(|key| key.as_str() != "name" && key.as_str() != SCHEMA_KEY)
    .collect();
  match tags.as_slice() {
    [] => diags.push(
      RawDiagnostic::new(
        "config::missing-command",
        format!("pipeline step has no command; expected one of {}", commands_list()),
      )
      .at(base),
    ),
    [tag] => {
      if !COMMAND_TAGS.contains(&tag.as_str()) {
        diags.push(
          RawDiagnostic::new("config::unknown-command", format!("unknown command `{tag}`"))
            .at(format!("{base}/{tag}"))
            .key_span()
            .help(suggestion_suffix(tag, &COMMAND_TAGS)),
        );
      }
    },
    _ => {
      let list = tags.iter().sorted().map(|tag| format!("`{tag}`")).join(", ");
      diags.push(
        RawDiagnostic::new(
          "config::multiple-commands",
          format!("pipeline step has more than one command ({list}); a step runs exactly one command"),
        )
        .at(base),
      );
    },
  }
}

pub fn interpolation_diagnostics(
  value: &Value,
  vars: &Map<String, Value>,
  step_names: &[String],
) -> Vec<RawDiagnostic> {
  let interp = Interpolator::default();
  let var_names: BTreeSet<&str> = vars.keys().map(String::as_str).collect();
  let step_index: BTreeMap<&str, usize> = step_names
    .iter()
    .enumerate()
    .map(|(position, name)| (name.as_str(), position))
    .collect();

  let mut diags = Vec::new();
  walk_leaves(value, "", &mut |pointer, leaf| {
    leaf_reference_diagnostics(&interp, pointer, leaf, &var_names, step_names, &step_index, &mut diags);
  });
  diags
}

fn leaf_reference_diagnostics(
  interp: &Interpolator,
  pointer: &str,
  leaf: &str,
  var_names: &BTreeSet<&str>,
  step_names: &[String],
  step_index: &BTreeMap<&str, usize>,
  diags: &mut Vec<RawDiagnostic>,
) {
  if !leaf.contains("{{") {
    return;
  }
  let scope = scope_of(pointer);

  for captures in step_ref().captures_iter(leaf) {
    let producer = &captures[1];
    match scope {
      Scope::Vars => diags.push(
        RawDiagnostic::new(
          "config::var-references-steps",
          "pipeline var references `steps`; vars may only use `vars` and `env`",
        )
        .at(pointer.to_owned()),
      ),
      Scope::Step(current) => match step_index.get(producer) {
        Some(&earlier) if earlier < current => {},
        Some(_) => diags.push(
          RawDiagnostic::new(
            "config::step-reference",
            format!(
              "step references `{producer}`, which is not an earlier step; steps may only reference steps before them"
            ),
          )
          .at(pointer.to_owned()),
        ),
        None => {
          let earlier: Vec<&str> = step_names.iter().take(current).map(String::as_str).collect();
          diags.push(
            RawDiagnostic::new(
              "config::step-reference",
              format!(
                "step references unknown step `{producer}`; {}",
                suggestion_suffix(producer, &earlier)
              ),
            )
            .at(pointer.to_owned()),
          );
        },
      },
      Scope::Other => diags.push(
        RawDiagnostic::new(
          "config::step-reference",
          "`steps` is not available here; only `vars` and `env` are",
        )
        .at(pointer.to_owned()),
      ),
    }
  }

  let stripped = step_ref().replace_all(leaf, "_").into_owned();
  match interp.references(&stripped) {
    Ok(references) => {
      for reference in references {
        let mut segments = reference.split('.');
        match segments.next() {
          Some("vars") => {
            if let Some(name) = segments.next() {
              if !var_names.contains(name) {
                let candidates: Vec<&str> = var_names.iter().copied().collect();
                diags.push(
                  RawDiagnostic::new(
                    "config::unknown-var",
                    format!("unknown variable `{name}`; {}", suggestion_suffix(name, &candidates)),
                  )
                  .at(pointer.to_owned()),
                );
              }
            }
          },
          Some("env" | "steps") => {},
          Some(other) => diags.push(
            RawDiagnostic::new(
              "config::unknown-namespace",
              format!("unknown namespace `{other}`; {}", suggestion_suffix(other, &NAMESPACES)),
            )
            .at(pointer.to_owned()),
          ),
          None => {},
        }
      }
    },
    Err(err) => {
      diags.push(RawDiagnostic::new("config::template", format!("invalid template: {err}")).at(pointer.to_owned()));
    },
  }
}

enum Scope {
  Vars,
  Step(usize),
  Other,
}

#[cfg_attr(
  dylint_lib = "treetime_lints",
  expect(
    result_defaulted,
    reason = "a path segment that is not an index places the pointer outside a step"
  )
)]
fn scope_of(pointer: &str) -> Scope {
  if let Some(rest) = pointer.strip_prefix("/steps/") {
    if let Some(index) = rest.split('/').next().and_then(|segment| segment.parse::<usize>().ok()) {
      return Scope::Step(index);
    }
  }
  if pointer == "/vars" || pointer.starts_with("/vars/") {
    return Scope::Vars;
  }
  Scope::Other
}

pub fn config_vars(value: &Value) -> Map<String, Value> {
  value
    .get("vars")
    .and_then(Value::as_object)
    .cloned()
    .unwrap_or_default()
}

pub fn config_step_names(value: &Value) -> Vec<String> {
  value
    .get("steps")
    .and_then(Value::as_array)
    .map(|steps| {
      steps
        .iter()
        .filter_map(|step| step.get("name").and_then(Value::as_str).map(str::to_owned))
        .collect()
    })
    .unwrap_or_default()
}

fn walk_leaves(value: &Value, pointer: &str, visit: &mut impl FnMut(&str, &str)) {
  match value {
    Value::String(text) => visit(pointer, text),
    Value::Array(items) => {
      for (position, item) in items.iter().enumerate() {
        walk_leaves(item, &format!("{pointer}/{position}"), visit);
      }
    },
    Value::Object(map) => {
      for (key, child) in map {
        walk_leaves(child, &format!("{pointer}/{}", escape_pointer(key)), visit);
      }
    },
    _ => {},
  }
}

fn string_options(options: &Value) -> Vec<String> {
  options
    .as_array()
    .map(|values| {
      values
        .iter()
        .filter_map(|value| value.as_str().map(str::to_owned))
        .collect()
    })
    .unwrap_or_default()
}

fn instance_string(value: &Value) -> String {
  match value {
    Value::String(text) => text.clone(),
    other => other.to_string(),
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use serde_json::json;

  fn codes(diags: &[RawDiagnostic]) -> Vec<String> {
    diags.iter().map(|diag| diag.code.clone()).sorted().collect()
  }

  fn find<'a>(diags: &'a [RawDiagnostic], code: &str) -> &'a RawDiagnostic {
    diags
      .iter()
      .find(|diag| diag.code == code)
      .expect("diagnostic with code present")
  }

  #[test]
  fn test_checks_structural_batches_all_shape_errors() {
    let value = json!({
      "bogus": 1,
      "steps": [
        { "name": "x", "timetree": {} },
        { "name": "x", "clock": {} },
        { "name": "y" },
        { "name": "z", "a": {}, "b": {} },
        { "timetree": {} },
        { "name": "w", "timtree": {} }
      ]
    });
    let diags = pipeline_structural_diagnostics(&value);
    assert_eq!(
      vec![
        "config::duplicate-step",
        "config::missing-command",
        "config::missing-name",
        "config::multiple-commands",
        "config::unknown-command",
        "config::unknown-field",
      ],
      codes(&diags)
    );
  }

  #[test]
  fn test_checks_structural_unknown_top_level_key_suggests() {
    let value = json!({ "step": [], "steps": [] });
    let diags = pipeline_structural_diagnostics(&value);
    let diag = find(&diags, "config::unknown-field");
    assert_eq!(Some("/step".to_owned()), diag.pointer);
    assert!(diag.use_key_span);
    assert_eq!(
      Some("did you mean `steps`? Valid values: `$schema`, `output_all`, `steps`, `vars`".to_owned()),
      diag.help
    );
  }

  #[test]
  fn test_checks_interpolation_batches_reference_errors() {
    let value = json!({
      "vars": { "data": "d", "slot": "{{ steps.tt.outputs.nwk }}" },
      "steps": [
        { "name": "tt", "timetree": { "tree": "{{ vars.dataa }}/t.nwk" } },
        { "name": "anc", "ancestral": { "tree": "{{ steps.zzz.outputs.nwk }}", "aln": "{{ foo.bar }}" } }
      ]
    });
    let vars = value["vars"].as_object().unwrap().clone();
    let step_names = vec!["tt".to_owned(), "anc".to_owned()];
    let diags = interpolation_diagnostics(&value, &vars, &step_names);
    assert_eq!(
      vec![
        "config::step-reference",
        "config::unknown-namespace",
        "config::unknown-var",
        "config::var-references-steps",
      ],
      codes(&diags)
    );
    assert_eq!(
      Some("/steps/0/timetree/tree".to_owned()),
      find(&diags, "config::unknown-var").pointer
    );
  }

  #[test]
  fn test_checks_interpolation_allows_backward_step_reference() {
    let value = json!({
      "steps": [
        { "name": "tt", "timetree": { "tree": "in.nwk" } },
        { "name": "anc", "ancestral": { "tree": "{{ steps.tt.outputs.nwk }}" } }
      ]
    });
    let diags = interpolation_diagnostics(&value, &Map::new(), &["tt".to_owned(), "anc".to_owned()]);
    assert!(diags.is_empty(), "a backward step reference must not be flagged");
  }

  #[test]
  fn test_checks_schema_reports_type_and_enum() {
    let value = json!({
      "steps": [
        { "name": "a", "timetree": { "clock_rate": "abc" } },
        { "name": "b", "ancestral": { "tree": "t.nwk", "method_anc": "bogus" } }
      ]
    });
    let diags = pipeline_schema_diagnostics(&value);
    let present = codes(&diags);
    assert!(
      present.contains(&"config::type".to_owned()),
      "type error expected, got {present:?}"
    );
    assert!(
      present.contains(&"config::enum".to_owned()),
      "enum error expected, got {present:?}"
    );
    assert_eq!(
      Some("/steps/1/ancestral/method_anc".to_owned()),
      find(&diags, "config::enum").pointer
    );
  }

  #[test]
  fn test_checks_schema_skips_template_in_typed_field() {
    let value = json!({
      "steps": [ { "name": "a", "timetree": { "clock_rate": "{{ vars.rate }}" } } ]
    });
    let diags = pipeline_schema_diagnostics(&value);
    assert!(
      diags.is_empty(),
      "template in a numeric field must be skipped, got {:?}",
      codes(&diags)
    );
  }
}
