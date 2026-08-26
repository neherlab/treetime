use crate::io::compression::remove_compression_ext;
use crate::io::file::{create_file_or_stdout, open_file_or_stdin};
use crate::io::fs::extension;
use crate::io::yaml::{yaml_read_file, yaml_write_file};
use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};
use serde_json::{Deserializer, de::Read};
use std::io::{Cursor, Write};
use std::path::Path;

pub fn json_read_file<T: for<'de> Deserialize<'de>, P: AsRef<Path>>(filepath: P) -> Result<T, Report> {
  let filepath = filepath.as_ref();
  json_read(open_file_or_stdin(&Some(filepath))?)
    .wrap_err_with(|| format!("When reading JSON file: '{}'", filepath.display()))
}

pub fn json_read_str<T: for<'de> Deserialize<'de>>(s: impl AsRef<str>) -> Result<T, Report> {
  json_read(Cursor::new(s.as_ref())).wrap_err("When reading JSON string")
}

pub fn json_read<T: for<'de> Deserialize<'de>>(reader: impl std::io::Read) -> Result<T, Report> {
  let mut de = Deserializer::from_reader(reader);
  deserialize_without_recursion_limit(&mut de).wrap_err("When reading JSON")
}

/// Mitigates recursion limit error when parsing large JSONs
/// See https://github.com/serde-rs/json/issues/334
fn deserialize_without_recursion_limit<'de, R: Read<'de>, T: Deserialize<'de>>(
  de: &mut Deserializer<R>,
) -> Result<T, Report> {
  de.disable_recursion_limit();
  let de = serde_stacker::Deserializer::new(de);
  let obj = T::deserialize(de).wrap_err("When parsing JSON")?;
  Ok(obj)
}

#[derive(Clone, Copy, Debug)]
pub struct JsonPretty(pub bool);

pub fn json_write_file<T: Serialize>(filepath: impl AsRef<Path>, obj: &T, pretty: JsonPretty) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  json_write(create_file_or_stdout(filepath)?, &obj, pretty)
    .wrap_err_with(|| format!("When writing JSON file: '{}'", filepath.display()))
}

pub fn json_write_str<T: Serialize>(obj: &T, pretty: JsonPretty) -> Result<String, Report> {
  if pretty.0 {
    serde_json::to_string_pretty(obj)
  } else {
    serde_json::to_string(obj)
  }
  .wrap_err("When writing JSON string")
}

pub fn json_write<W: Write, T: Serialize>(writer: W, obj: &T, pretty: JsonPretty) -> Result<(), Report> {
  if pretty.0 {
    serde_json::to_writer_pretty(writer, &obj)
  } else {
    serde_json::to_writer(writer, &obj)
  }
  .wrap_err("When writing JSON")
}

/// Writes JSON or YAML file depending on file extension
pub fn json_or_yaml_write_file<T: Serialize>(filepath: impl AsRef<Path>, obj: &T) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  match serialization_format_from_path(filepath) {
    SerializationFormat::Json => json_write_file(filepath, &obj, JsonPretty(true)),
    SerializationFormat::Yaml => yaml_write_file(filepath, &obj),
  }
}

/// Reads JSON or YAML file depending on file extension.
///
/// The format is chosen from the extension after any compression suffix is stripped (`.yaml`/`.yml`
/// select YAML, everything else JSON), mirroring `json_or_yaml_write_file`. Both readers decompress
/// transparently and accept `-` for stdin.
pub fn json_or_yaml_read_file<T: for<'de> Deserialize<'de>>(filepath: impl AsRef<Path>) -> Result<T, Report> {
  let filepath = filepath.as_ref();
  match serialization_format_from_path(filepath) {
    SerializationFormat::Json => json_read_file(filepath),
    SerializationFormat::Yaml => yaml_read_file(filepath),
  }
}

pub(crate) fn serialization_format_from_path(filepath: impl AsRef<Path>) -> SerializationFormat {
  let filepath = remove_compression_ext(filepath);
  match extension(filepath).map(|ext| ext.to_lowercase()) {
    Some(ext) if matches!(ext.as_str(), "yaml" | "yml") => SerializationFormat::Yaml,
    _ => SerializationFormat::Json,
  }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum SerializationFormat {
  Json,
  Yaml,
}

/// Check whether a serde value serializes to null. This is useful to skip a generic struct field even if we don't
/// know the exact type. Usage: add attribute `#[serde(skip_serializing_if = "is_json_value_null")]` to a struct field
/// you want to skip.
pub fn is_json_value_null<T: Serialize>(t: &T) -> bool {
  match serde_json::to_value(t) {
    Ok(v) => v.is_null(),
    Err(e) => {
      log::warn!("JSON serialization failed during null check, treating as null: {e}");
      true
    },
  }
}
