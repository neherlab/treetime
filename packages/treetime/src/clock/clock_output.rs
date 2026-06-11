use crate::clock::clock_model::ClockModel;
use eyre::Report;
use std::path::Path;
use treetime_utils::io::json::{JsonPretty, json_write_file};

/// Write clock model parameters to a JSON file at the given path.
pub fn write_clock_model(clock_model: &ClockModel, path: &Path) -> Result<(), Report> {
  json_write_file(path, clock_model, JsonPretty(true))
}
