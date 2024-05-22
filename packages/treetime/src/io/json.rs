use crate::io::file::create_file;
use eyre::{Report, WrapErr};
use serde::Serialize;
use std::io::Write;
use std::path::Path;

#[derive(Clone, Copy, Debug)]
pub struct JsonPretty(pub bool);

pub fn json_stringify<T: Serialize>(obj: &T, pretty: JsonPretty) -> Result<String, Report> {
  if pretty.0 {
    serde_json::to_string_pretty(obj)
  } else {
    serde_json::to_string(obj)
  }
  .wrap_err("When converting an entry to JSON string")
}

pub fn json_write_file<T: Serialize>(obj: &T, filepath: impl AsRef<Path>, pretty: JsonPretty) -> Result<(), Report> {
  json_write_writer(&obj, &mut create_file(filepath)?, pretty)
}

pub fn json_write_writer<W: Write, T: Serialize>(obj: &T, writer: W, pretty: JsonPretty) -> Result<(), Report> {
  if pretty.0 {
    serde_json::to_writer_pretty(writer, &obj)
  } else {
    serde_json::to_writer(writer, &obj)
  }
  .wrap_err("When writing JSON")
}
