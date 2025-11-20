use crate::commands::clock::clock_model::ClockModel;
use eyre::Report;
use std::path::Path;
use treetime_io::json::{JsonPretty, json_write_file};

/// Write clock model parameters to file.
///
/// Core: Documents inferred molecular clock rate.
/// Why: Clock rate is key parameter for interpreting divergence times.
/// How: JSON file with rate, intercept, R², confidence intervals.
pub fn write_clock_model(clock_model: &ClockModel, out_base: &Path) -> Result<(), Report> {
  let json_path = out_base.with_extension("json");
  json_write_file(json_path, clock_model, JsonPretty(true))
}
