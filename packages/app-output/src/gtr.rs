use eyre::Report;
use std::path::Path;
use treetime::gtr::get_gtr::GtrOutput;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn write_gtr_json(output: &GtrOutput, path: impl AsRef<Path>) -> Result<(), Report> {
  json_write_file(path, output, JsonPretty(true))
}
