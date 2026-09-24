use crate::config::source::RawDiagnostic;
use crate::config::suggest::suggestion_suffix;
use jsonschema::error::ValidationErrorKind;
use schemars::Schema;
use serde_json::Value;

pub fn schema_diagnostics(value: &Value, schema: &Schema, skip_templates: bool) -> Vec<RawDiagnostic> {
  schema_diagnostics_prefixed(value, schema, skip_templates, "")
}

pub fn schema_diagnostics_prefixed(
  value: &Value,
  schema: &Schema,
  skip_templates: bool,
  base: &str,
) -> Vec<RawDiagnostic> {
  let schema_value = match serde_json::to_value(schema) {
    Ok(schema_value) => schema_value,
    Err(err) => {
      return vec![RawDiagnostic::builder("config::internal", format!("could not build schema: {err}")).build()];
    },
  };
  let validator = match jsonschema::validator_for(&schema_value) {
    Ok(validator) => validator,
    Err(err) => return vec![RawDiagnostic::builder("config::internal", format!("invalid schema: {err}")).build()],
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
          RawDiagnostic::builder("config::enum", format!("`{bad}` is not a valid value"))
            .at(pointer)
            .help(suggestion_suffix(&bad, &candidate_refs))
            .build(),
        );
      },
      ValidationErrorKind::Type { .. } => {
        diags.push(
          RawDiagnostic::builder("config::type", error.to_string())
            .at(pointer)
            .build(),
        );
      },
      ValidationErrorKind::Required { property } => {
        let property = property.as_str().unwrap_or_default();
        diags.push(
          RawDiagnostic::builder("config::required", format!("missing required field `{property}`"))
            .at(pointer)
            .build(),
        );
      },
      ValidationErrorKind::AdditionalProperties { unexpected } => {
        for key in unexpected {
          diags.push(
            RawDiagnostic::builder("config::unknown-field", format!("unknown field `{key}`"))
              .at(format!("{pointer}/{key}"))
              .key_span(true)
              .build(),
          );
        }
      },
      _ => {
        diags.push(
          RawDiagnostic::builder("config::schema", error.to_string())
            .at(pointer)
            .build(),
        );
      },
    }
  }
  diags
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
