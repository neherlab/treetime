use crate::cli::pipeline::types::PipelineStepCommand;
use serde_json::Value;

/// Input fields a step reads, as (label, key path into the command's args object).
///
/// One source of truth for both the dry-run plan and the plan-time safety checks, so the two never
/// disagree about what a step reads.
pub const INPUT_FIELDS: [(&str, &[&str]); 5] = [
  ("tree", &["tree"]),
  ("alignment", &["alignment"]),
  ("metadata", &["metadata"]),
  ("weights", &["weights"]),
  ("vcf-reference", &["vcf_reference"]),
];

/// Labeled input paths a step reads, for the dry-run plan.
pub fn labeled_input_paths(command: &PipelineStepCommand) -> Vec<(&'static str, String)> {
  let args = command.args_value();
  INPUT_FIELDS
    .iter()
    .flat_map(|(label, keys)| {
      lookup_paths(&args, keys)
        .into_iter()
        .map(move |path| (*label, path.to_owned()))
    })
    .collect()
}

/// Every input path a step reads, unlabeled, for the safety checks.
pub fn input_paths(command: &PipelineStepCommand) -> Vec<String> {
  labeled_input_paths(command).into_iter().map(|(_, path)| path).collect()
}

/// Follow a key path into the args object and return the string paths it holds.
///
/// A leaf may be a single string (`tree`, `metadata`) or a list of strings (`alignment` accepts
/// several files); both are flattened to individual paths.
fn lookup_paths<'a>(args: &'a Value, keys: &[&str]) -> Vec<&'a str> {
  let mut current = args;
  for key in keys {
    let Some(next) = current.get(key) else {
      return Vec::new();
    };
    current = next;
  }
  match current {
    Value::String(path) => vec![path.as_str()],
    Value::Array(items) => items.iter().filter_map(Value::as_str).collect(),
    _ => Vec::new(),
  }
}
