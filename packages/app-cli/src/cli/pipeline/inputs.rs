use crate::cli::pipeline::types::PipelineStepCommand;
use eyre::Report;
use serde_json::Value;

const INPUT_FIELDS: [(&str, &[&str]); 5] = [
  ("tree", &["tree"]),
  ("alignment", &["alignment"]),
  ("metadata", &["metadata"]),
  ("weights", &["weights"]),
  ("vcf-reference", &["vcf_reference"]),
];

pub fn input_paths(command: &PipelineStepCommand) -> Result<Vec<String>, Report> {
  Ok(
    labeled_input_paths(command)?
      .into_iter()
      .map(|(_, path)| path)
      .collect(),
  )
}

pub fn labeled_input_paths(command: &PipelineStepCommand) -> Result<Vec<(&'static str, String)>, Report> {
  let args = command.args_value()?;
  Ok(
    INPUT_FIELDS
      .iter()
      .flat_map(|(label, keys)| {
        lookup_paths(&args, keys)
          .into_iter()
          .map(move |path| (*label, path.to_owned()))
      })
      .collect(),
  )
}

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
