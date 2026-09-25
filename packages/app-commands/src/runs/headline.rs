use crate::command::{AppCommand, OutputFile};
use crate::json_float::JsonFloat;
use app_output::output_plan::OutputSelection;
use eyre::{Report, WrapErr};
use serde_json::Value;
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

pub const HEADLINE_CLOCK_RATE: &str = "clock_rate";
pub const HEADLINE_R_SQUARED: &str = "r_squared";
pub const HEADLINE_ROOT_DATE: &str = "root_date";

pub fn run_headline(
  command: AppCommand,
  out_dir: &Path,
  output_files: &[OutputFile],
) -> Result<BTreeMap<String, JsonFloat>, Report> {
  let mut headline = BTreeMap::new();
  for file in output_files {
    let path = out_dir.join(&file.path);
    match file.kind {
      OutputSelection::ClockModel => {
        let model = read_json(&path)?;
        if let Some(rate) = model.get("clock_rate").and_then(Value::as_f64) {
          headline.insert(HEADLINE_CLOCK_RATE.to_owned(), JsonFloat(rate));
        }
        if let Some(r) = model.pointer("/stats/estimated/r_val").and_then(Value::as_f64) {
          headline.insert(HEADLINE_R_SQUARED.to_owned(), JsonFloat(r * r));
        }
      },
      OutputSelection::Auspice if command == AppCommand::Timetree => {
        let tree = read_json(&path)?;
        if let Some(date) = tree.pointer("/tree/node_attrs/num_date/value").and_then(Value::as_f64) {
          headline.insert(HEADLINE_ROOT_DATE.to_owned(), JsonFloat(date));
        }
      },
      _ => {},
    }
  }
  Ok(headline)
}

fn read_json(path: &Path) -> Result<Value, Report> {
  let text = fs::read_to_string(path).wrap_err_with(|| format!("When reading '{}'", path.display()))?;
  serde_json::from_str(&text).wrap_err_with(|| format!("When parsing '{}'", path.display()))
}
