use eyre::{Report, WrapErr};
use serde::Serialize;
use std::io::Write;

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

pub fn json_write_impl<W: Write, T: Serialize>(writer: W, obj: &T, pretty: JsonPretty) -> Result<(), Report> {
  if pretty.0 {
    serde_json::to_writer_pretty(writer, &obj)
  } else {
    serde_json::to_writer(writer, &obj)
  }
  .wrap_err("When writing JSON")
}
