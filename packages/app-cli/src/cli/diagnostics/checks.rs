use crate::cli::diagnostics::source::RawDiagnostic;
use crate::cli::pipeline::interpolate::{Interpolator, NAMESPACES};
use crate::cli::pipeline::resolve::{STEP_REF, TOP_LEVEL_KEYS};
use crate::cli::pipeline::suggest::suggestion_suffix;
use crate::cli::pipeline::types::{COMMAND_TAGS, SCHEMA_KEY};
use crate::cli::schema::command_schema_for;
use itertools::Itertools;
use jsonschema::error::ValidationErrorKind;
use schemars::Schema;
use serde_json::{Map, Value};
use std::collections::{BTreeMap, BTreeSet};

/// Collect precise schema violations for every step's command payload.
///
/// Each step's command object is validated against that command's own strict schema rather than the
/// whole-pipeline schema, so a bad `enum`/`type` value keeps its exact location and message instead of
/// collapsing into the step's `oneOf` branch failure. Instance paths are prefixed with the step path
/// so a caret lands on the offending leaf. A `{{ ... }}` template in a typed field is skipped: it is an
/// interpolation slot resolved later, not a literal that must match the field type.
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

/// Collect every JSON-Schema violation of `value` against `schema` as a diagnostic.
///
/// `jsonschema` reports all structural problems in one pass, which is what makes the batch real
/// (type, enum, required, and unknown-field errors together). With `skip_templates`, a leaf that
/// holds a `{{ ... }}` template is left to the interpolation pass: it is an interpolation slot, not a
/// literal that must match the field's type.
pub fn schema_diagnostics(value: &Value, schema: &Schema, skip_templates: bool) -> Vec<RawDiagnostic> {
  schema_diagnostics_prefixed(value, schema, skip_templates, "")
}

/// Validate `value` against `schema`, prefixing every reported pointer with `base`.
///
/// `base` is empty for a whole-document validation and `/steps/<i>/<tag>` when validating a step's
/// command payload in isolation, which keeps carets pointing at the original file location.
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

/// Collect the pipeline's shape errors (keys, step names, command tags) as diagnostics.
///
/// Mirrors the fail-fast checks in the loader, but gathers every problem so the whole document is
/// reported at once. Shape errors are checked before schema and interpolation because a malformed
/// skeleton (missing `steps`, a non-mapping step) makes those later passes noisy.
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

/// Collect one step's shape errors: name presence and type, duplicate names, and the command tag.
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

/// Collect static interpolation-reference errors for every string leaf in the document.
///
/// This is a static check: it validates that each `{{ ... }}` reference names a known namespace, a
/// declared `vars` entry, or an earlier step, and enforces backward-only step references. Output
/// selections and `env` values are resolved later, so they are not checked here.
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

/// Validate the references in one string leaf, given the leaf's JSON pointer for scope and carets.
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

  for captures in STEP_REF.captures_iter(leaf) {
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

  let stripped = STEP_REF.replace_all(leaf, "_").into_owned();
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
      diags.push(RawDiagnostic::new("config::template", format!("invalid template: {err}")).at(pointer.to_owned()))
    },
  }
}

/// The interpolation scope a leaf sits in, which governs whether `steps` references are allowed.
enum Scope {
  Vars,
  Step(usize),
  Other,
}

/// Determine a leaf's interpolation scope from its JSON pointer.
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

/// The `vars` mapping of a config value, or an empty map when absent or mistyped.
pub fn config_vars(value: &Value) -> Map<String, Value> {
  value
    .get("vars")
    .and_then(Value::as_object)
    .cloned()
    .unwrap_or_default()
}

/// The step names of a config value, in list order, skipping steps without a string name.
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

/// Visit every string leaf of a value, passing its JSON pointer and text to `visit`.
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
        walk_leaves(child, &format!("{pointer}/{}", escape_segment(key)), visit);
      }
    },
    _ => {},
  }
}

/// Escape a mapping key for use as a JSON-pointer segment (RFC 6901).
fn escape_segment(segment: &str) -> String {
  segment.replace('~', "~0").replace('/', "~1")
}

/// Human-readable list of the accepted command tags, in declared order.
fn commands_list() -> String {
  COMMAND_TAGS.iter().map(|tag| format!("`{tag}`")).join(", ")
}

/// Extract the string entries of a JSON array (used for enum option lists).
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

/// Render a JSON instance as a short string for an error message.
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

  // The structural pass gathers every shape problem in one run rather than stopping at the first.
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

  // An unknown top-level key is reported at its own pointer with a did-you-mean over the valid keys.
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

  // The interpolation pass batches every static reference error, each pointed at its own leaf.
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

  // A backward step reference (earlier step) is accepted; a forward one is not.
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

  // Per-step schema validation reports precise type and enum errors on literal command values, each
  // located at the offending leaf rather than collapsed into the step's command `oneOf`.
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

  // A whole-value template in a typed field is an interpolation slot, not a schema type error.
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
