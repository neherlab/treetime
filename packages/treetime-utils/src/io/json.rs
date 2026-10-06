use crate::io::file::{read_file_with, write_file_with};
use eyre::{Report, WrapErr};
use serde::Serialize;
use serde::de::DeserializeOwned;
use serde_json::Deserializer;
use std::io::{Read, Write};
use std::path::Path;

pub fn json_read_file<T: DeserializeOwned>(filepath: impl AsRef<Path>) -> Result<T, Report> {
  read_file_with(filepath, json_read)
}

pub fn json_read_str<T: DeserializeOwned>(s: impl AsRef<str>) -> Result<T, Report> {
  json_read(s.as_ref().as_bytes())
}

pub fn json_read<T: DeserializeOwned>(reader: impl Read) -> Result<T, Report> {
  let mut de = Deserializer::from_reader(reader);
  de.disable_recursion_limit();
  let value = T::deserialize(serde_stacker::Deserializer::new(&mut de)).wrap_err("When parsing JSON")?;
  de.end().wrap_err("When parsing JSON")?;
  Ok(value)
}

pub fn json_read_slice<T: DeserializeOwned>(bytes: &[u8]) -> Result<T, Report> {
  let mut de = Deserializer::from_slice(bytes);
  de.disable_recursion_limit();
  let value = T::deserialize(serde_stacker::Deserializer::new(&mut de)).wrap_err("When parsing JSON")?;
  de.end().wrap_err("When parsing JSON")?;
  Ok(value)
}

pub fn json_write_file<T: Serialize>(filepath: impl AsRef<Path>, obj: &T, pretty: JsonPretty) -> Result<(), Report> {
  write_file_with(filepath, |writer| {
    json_write(&mut *writer, obj, pretty)?;
    writeln!(writer)?;
    Ok(())
  })
}

pub fn json_write_str<T: Serialize>(obj: &T, pretty: JsonPretty) -> Result<String, Report> {
  let mut buf = Vec::new();
  json_write(&mut buf, obj, pretty)?;
  Ok(String::from_utf8(buf)?)
}

pub fn json_write<T: Serialize>(writer: impl Write, obj: &T, pretty: JsonPretty) -> Result<(), Report> {
  if pretty.0 {
    serde_json::to_writer_pretty(writer, obj)
  } else {
    serde_json::to_writer(writer, obj)
  }
  .wrap_err("When writing JSON")
}

#[derive(Clone, Copy, Debug)]
pub struct JsonPretty(pub bool);

pub fn is_json_value_null<T: Serialize>(t: &T) -> bool {
  match serde_json::to_value(t) {
    Ok(v) => v.is_null(),
    Err(e) => {
      log::warn!("JSON serialization failed during null check, treating as null: {e}");
      true
    },
  }
}
