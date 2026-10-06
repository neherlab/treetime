use crate::io::file::{read_file_with, write_file_with};
use crate::make_report;
use deser::Serialize;
use deser::adapters::As;
use deser::de::DeserializeOwned;
use deser_json::{Deserializer, DeserializerConfig, Indent, SerializerConfig};
use deser_path::PathLayer;
use deser_serde::Serde;
use eyre::{Report, WrapErr};
use serde_json::Value;
use std::io::{Read, Write};
use std::path::Path;

const COMPACT: SerializerConfig = SerializerConfig::new();

const PRETTY: SerializerConfig = SerializerConfig::new().pretty(Indent::Spaces(2));

pub fn json_read_file<T: DeserializeOwned>(filepath: impl AsRef<Path>) -> Result<T, Report> {
  read_file_with(filepath, json_read)
}

pub fn json_read_str<T: DeserializeOwned>(s: impl AsRef<str>) -> Result<T, Report> {
  Deserializer::from_str(s.as_ref())
    .deserialize_with(|driver| driver.push_layer(PathLayer::new()))
    .wrap_err("When parsing JSON")
}

pub fn json_read<T: DeserializeOwned>(reader: impl Read) -> Result<T, Report> {
  let mut reader = DeserializerConfig::new().reader(reader);
  let value = reader
    .read_with(|driver| driver.push_layer(PathLayer::new()))
    .wrap_err("When parsing JSON")?
    .ok_or_else(|| make_report!("When parsing JSON: the input is empty"))?;
  reader.end().wrap_err("When parsing JSON")?;
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
  config(pretty).to_string(obj).wrap_err("When writing JSON")
}

pub fn json_write<T: Serialize>(writer: impl Write, obj: &T, pretty: JsonPretty) -> Result<(), Report> {
  config(pretty).to_writer(writer, obj).wrap_err("When writing JSON")
}

pub fn json_value_read_str(s: impl AsRef<str>) -> Result<Value, Report> {
  json_read_str::<As<Value, Serde>>(s).map(As::into_inner)
}

pub fn json_value_read_file(filepath: impl AsRef<Path>) -> Result<Value, Report> {
  json_read_file::<As<Value, Serde>>(filepath).map(As::into_inner)
}

pub fn to_json_value<T: Serialize>(obj: &T) -> Result<Value, Report> {
  let value = deser_value::to_value(obj).wrap_err("When converting to a JSON value")?;
  let value: As<Value, Serde> = deser_value::from_value(&value).wrap_err("When converting to a JSON value")?;
  Ok(value.into_inner())
}

pub fn from_json_value<T: DeserializeOwned>(value: &Value) -> Result<T, Report> {
  let value = deser_value::to_value(&As::<_, Serde>::new(value)).wrap_err("When converting a JSON value")?;
  deser_value::from_value(&value).wrap_err("When converting a JSON value")
}

#[derive(Clone, Copy, Debug)]
pub struct JsonPretty(pub bool);

const fn config(pretty: JsonPretty) -> &'static SerializerConfig {
  if pretty.0 { &PRETTY } else { &COMPACT }
}
