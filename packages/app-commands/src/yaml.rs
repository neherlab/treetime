use deser::adapters::As;
use deser::de::{DeserializeOwned, Layer, LayerEvent, Next};
use deser::{Atom, Error, Event};
use deser_path::PathLayer;
use deser_serde::Serde;
use deser_value::Kind;
use deser_yaml::{Deserializer, DeserializerConfig};
use eyre::Report;
use itertools::Itertools;
use saphyr::{Mapping, Scalar, ScalarStyle, Yaml, YamlEmitter};
use serde_json::Value;
use std::borrow::Cow;
use std::sync::Arc;
use treetime_utils::io::json::{JsonPretty, json_write_str};
use treetime_utils::{make_error, make_report};

const DOCUMENT_START: &str = "---\n";

const YAML: DeserializerConfig = DeserializerConfig::new().track_locations(true);

pub fn yaml_read_str<T: DeserializeOwned>(text: &str) -> Result<T, Report> {
  Deserializer::from_str_with_config(text, &YAML)
    .deserialize_with(|driver| {
      driver.push_layer(PathLayer::new());
      driver.push_layer(PlainBooleanWords { text: Arc::from(text) });
    })
    .map_err(Report::new)
}

pub fn yaml_value_read_str(text: &str) -> Result<Value, Report> {
  let value: deser_value::Value = yaml_read_str(text)?;
  ensure_finite_numbers(&value)?;
  let value: As<Value, Serde> = deser_value::from_value(&value).map_err(Report::new)?;
  Ok(value.into_inner())
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

fn ensure_finite_numbers(value: &deser_value::Value) -> Result<(), Report> {
  match value.kind() {
    Kind::Seq(items) => items.iter().try_for_each(ensure_finite_numbers),
    Kind::Map(entries) => entries.iter().try_for_each(|(_, item)| ensure_finite_numbers(item)),
    kind => match kind.as_f64() {
      Some(number) if !number.is_finite() => {
        let location = value.span().map(|span| span.start()).map_or_else(String::new, |start| {
          format!(" at line {} column {}", start.line, start.column)
        });
        make_error!("the number {number} has no JSON representation{location}")
      },
      _ => Ok(()),
    },
  }
}

struct PlainBooleanWords {
  text: Arc<str>,
}

impl Layer for PlainBooleanWords {
  fn event<'de>(&mut self, event: LayerEvent<'_, 'de>, next: &mut Next<'_, 'de>) -> Result<(), Error> {
    let word = match event.event() {
      Event::Atom(Atom::Str(text)) if !next.state().is_map_key() => boolean_word(text),
      _ => None,
    };
    let plain = next
      .state()
      .input_range()
      .and_then(|range| self.text.get(range))
      .is_some_and(|source| !source.starts_with(['"', '\'', '|', '>']));
    match word {
      Some(word) if plain => next.emit(LayerEvent::new(Event::Atom(Atom::Bool(word)))),
      _ => next.emit(event),
    }
  }
}

fn boolean_word(text: &str) -> Option<bool> {
  const TRUE: [&str; 3] = ["y", "yes", "on"];
  const FALSE: [&str; 3] = ["n", "no", "off"];
  if TRUE.iter().any(|word| text.eq_ignore_ascii_case(word)) {
    Some(true)
  } else if FALSE.iter().any(|word| text.eq_ignore_ascii_case(word)) {
    Some(false)
  } else {
    None
  }
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
