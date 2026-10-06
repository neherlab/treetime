use deser::adapters::As;
use deser_serde::Serde;
use eyre::Report;
use itertools::Itertools;
use saphyr::{Mapping, Scalar, ScalarStyle, Yaml, YamlEmitter};
use serde::de::DeserializeOwned;
use serde_json::Value;
use serde_saphyr::DuplicateKeyPolicy;
use std::borrow::Cow;
use treetime_utils::io::json::{JsonPretty, json_write_str};
use treetime_utils::make_report;

const DOCUMENT_START: &str = "---\n";

pub fn yaml_read_str<T: DeserializeOwned>(text: &str) -> Result<T, Report> {
  let options = serde_saphyr::options! {
    duplicate_keys: DuplicateKeyPolicy::Error,
    reject_non_finite_typeless_float: true,
  };
  serde_saphyr::from_str_with_options(text, options).map_err(Report::new)
}

pub fn yaml_text(key: &str, value: &Value) -> Result<String, Report> {
  let mut entry = Mapping::new();
  entry.insert(
    Yaml::Value(Scalar::String(Cow::Owned(key.to_owned()))),
    yaml_node(value)?,
  );
  emit(&Yaml::Mapping(entry)).map_err(|err| make_report!("could not write `{key}` as YAML: {err}"))
}

pub fn yaml_document(value: &Value) -> Result<String, Report> {
  let text = emit(&yaml_node(value)?).map_err(|err| make_report!("could not write a YAML document: {err}"))?;
  Ok(format!("{text}\n"))
}

fn emit(node: &Yaml<'_>) -> Result<String, Report> {
  let mut text = String::new();
  YamlEmitter::new(&mut text).dump(node).map_err(Report::new)?;
  Ok(text.strip_prefix(DOCUMENT_START).unwrap_or(&text).to_owned())
}

fn yaml_node(value: &Value) -> Result<Yaml<'static>, Report> {
  let plain = |text: String| Yaml::Representation(Cow::Owned(text), ScalarStyle::Plain, None);
  Ok(match value {
    Value::Null => plain("null".to_owned()),
    Value::Bool(_) | Value::Number(_) => plain(value.to_string()),
    Value::String(_) => {
      let quoted = json_write_str(&As::<_, Serde>::new(value), JsonPretty(false))?;
      let escaped = quoted
        .strip_prefix('"')
        .and_then(|text| text.strip_suffix('"'))
        .ok_or_else(|| make_report!("a JSON string must be quoted: {quoted}"))?;
      Yaml::Representation(Cow::Owned(escaped.to_owned()), ScalarStyle::DoubleQuoted, None)
    },
    Value::Array(items) => Yaml::Sequence(items.iter().map(yaml_node).try_collect()?),
    Value::Object(entries) => {
      let mut mapping = Mapping::new();
      for (key, item) in entries {
        mapping.insert(Yaml::Value(Scalar::String(Cow::Owned(key.clone()))), yaml_node(item)?);
      }
      Yaml::Mapping(mapping)
    },
  })
}
