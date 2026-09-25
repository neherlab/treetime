use crate::command::AppCommand;
use crate::config::catalog::{CommandSettings, SettingRole, SettingSpec, command_settings};
use crate::config::settings::setting_ref;
use crate::run_config::RUN_CONFIG_OUTPUT_DIR;
use app_datasets::schema_directive;
use clap::{Arg, Command};
use eyre::Report;
use itertools::Itertools;
use saphyr::{Mapping, Scalar, ScalarStyle, Yaml, YamlEmitter};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::borrow::Cow;
use treetime_utils::io::json::{JsonPretty, json_write_str};
use treetime_utils::{make_error, make_report};

const OUTPUT_ALL_KEY: &str = "output_all";

const OUTPUT_ALL_FLAG: &str = "--output-all";

const DOCUMENT_START: &str = "---\n";

pub fn config_code(command: AppCommand, config: &Map<String, Value>) -> Result<ConfigCode, Report> {
  let settings = command_settings(command)?;
  let cli = command.cli_command();
  let inputs = settings
    .settings
    .iter()
    .filter(|spec| matches!(spec.role, SettingRole::Input | SettingRole::InputTemplate))
    .map(|spec| (spec, setting_value(config, spec)))
    .filter(|(_, value)| has_path(value))
    .collect_vec();
  let changed = changed_settings(&settings, config);
  let output_dir = config
    .get(OUTPUT_ALL_KEY)
    .and_then(Value::as_str)
    .unwrap_or(RUN_CONFIG_OUTPUT_DIR);

  let command_line = command_line(command, &cli, &inputs, &changed, output_dir)?;
  let yaml = yaml(command, &inputs, &changed, output_dir)?;
  Ok(ConfigCode {
    command_line_text: command_line_text(&command_line),
    yaml_text: format!("{}\n", yaml.iter().map(|line| &line.text).join("\n")),
    command_line,
    yaml,
  })
}

pub fn setting_tokens(arg: &Arg, flag: &str, value: &Value) -> Option<Vec<String>> {
  if !arg.get_action().takes_values() {
    return (value == &Value::Bool(true)).then(|| vec![flag.to_owned()]);
  }
  match value {
    Value::Null => None,
    Value::Array(items) => list_tokens(arg, flag, items),
    _ => Some(vec![flag.to_owned(), cli_value(arg, value)]),
  }
}

pub fn yaml_text(key: &str, value: &Value) -> Result<String, Report> {
  let mut entry = Mapping::new();
  entry.insert(
    Yaml::Value(Scalar::String(Cow::Owned(key.to_owned()))),
    yaml_node(value)?,
  );
  let mut text = String::new();
  YamlEmitter::new(&mut text)
    .dump(&Yaml::Mapping(entry))
    .map_err(|err| make_report!("could not write `{key}` as YAML: {err}"))?;
  Ok(text.strip_prefix(DOCUMENT_START).unwrap_or(&text).to_owned())
}

/// A command line and a YAML config that reproduce a configuration.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct ConfigCode {
  /// The command line, one flag per line.
  pub command_line: Vec<CodeLine>,
  /// The command line as a shell reads it, with a line continuation after every flag but the last.
  pub command_line_text: String,
  /// The YAML config, line by line.
  pub yaml: Vec<CodeLine>,
  /// The YAML config as a file holds it.
  pub yaml_text: String,
}

/// One line of a command line or a YAML config.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct CodeLine {
  /// Text of the line.
  pub text: String,
  /// What the line sets.
  pub kind: CodeLineKind,
}

/// What a line of a command line or a YAML config sets.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum CodeLineKind {
  /// The command.
  Command,
  /// An input file.
  Input,
  /// A setting that differs from its default.
  Changed,
  /// The output directory.
  Output,
  /// A comment.
  Comment,
}

fn changed_settings<'a>(settings: &'a CommandSettings, config: &Map<String, Value>) -> Vec<(&'a SettingSpec, Value)> {
  settings
    .settings
    .iter()
    .filter(|spec| spec.role == SettingRole::Setting)
    .map(|spec| (spec, setting_value(config, spec)))
    .filter(|(spec, value)| *value != spec.default_value)
    .collect()
}

fn setting_value(config: &Map<String, Value>, spec: &SettingSpec) -> Value {
  setting_ref(config, &spec.path)
    .cloned()
    .unwrap_or_else(|| spec.default_value.clone())
}

fn has_path(value: &Value) -> bool {
  match value {
    Value::String(path) => !path.is_empty(),
    Value::Array(paths) => paths.iter().any(has_path),
    _ => false,
  }
}

fn command_line(
  command: AppCommand,
  cli: &Command,
  inputs: &[(&SettingSpec, Value)],
  changed: &[(&SettingSpec, Value)],
  output_dir: &str,
) -> Result<Vec<CodeLine>, Report> {
  let mut lines = vec![CodeLine {
    text: format!("treetime {command}"),
    kind: CodeLineKind::Command,
  }];
  for (settings, kind) in [(inputs, CodeLineKind::Input), (changed, CodeLineKind::Changed)] {
    for (spec, value) in settings {
      let arg = cli
        .get_arguments()
        .find(|arg| arg.get_long() == spec.flag.strip_prefix("--"))
        .ok_or_else(|| make_report!("setting `{}` has no command-line flag", spec.key))?;
      lines.push(match setting_tokens(arg, &spec.flag, value) {
        Some(tokens) => CodeLine {
          text: quote_tokens(&tokens)?,
          kind,
        },
        None => CodeLine {
          text: format!(
            "# {} = {} has no command-line form; use the YAML config",
            spec.key,
            json_write_str(value, JsonPretty(false))?
          ),
          kind: CodeLineKind::Comment,
        },
      });
    }
  }
  lines.push(CodeLine {
    text: quote_tokens(&[OUTPUT_ALL_FLAG.to_owned(), output_dir.to_owned()])?,
    kind: CodeLineKind::Output,
  });
  Ok(lines)
}

fn command_line_text(lines: &[CodeLine]) -> String {
  lines
    .iter()
    .filter(|line| line.kind != CodeLineKind::Comment)
    .enumerate()
    .map(|(index, line)| {
      if index == 0 {
        line.text.clone()
      } else {
        format!("  {}", line.text)
      }
    })
    .join(" \\\n")
}

fn yaml(
  command: AppCommand,
  inputs: &[(&SettingSpec, Value)],
  changed: &[(&SettingSpec, Value)],
  output_dir: &str,
) -> Result<Vec<CodeLine>, Report> {
  let mut lines = vec![
    CodeLine {
      text: schema_directive(command.into()),
      kind: CodeLineKind::Comment,
    },
    CodeLine {
      text: format!("# treetime {command} --config run.yaml"),
      kind: CodeLineKind::Comment,
    },
  ];
  for (spec, value) in inputs {
    lines.extend(yaml_entry(&spec.key, value, CodeLineKind::Input)?);
  }
  let mut nested = Map::new();
  for (spec, value) in changed {
    insert_at(&mut nested, &spec.path, value.clone())?;
  }
  for (key, value) in &nested {
    lines.extend(yaml_entry(key, value, CodeLineKind::Changed)?);
  }
  lines.extend(yaml_entry(
    OUTPUT_ALL_KEY,
    &Value::String(output_dir.to_owned()),
    CodeLineKind::Output,
  )?);
  Ok(lines)
}

fn yaml_entry(key: &str, value: &Value, kind: CodeLineKind) -> Result<Vec<CodeLine>, Report> {
  let text = yaml_text(key, value)?;
  Ok(
    text
      .lines()
      .map(|line| CodeLine {
        text: line.to_owned(),
        kind,
      })
      .collect(),
  )
}

fn insert_at(map: &mut Map<String, Value>, path: &[String], value: Value) -> Result<(), Report> {
  let Some((last, parents)) = path.split_last() else {
    return make_error!("a setting must have a key");
  };
  let mut current = map;
  for key in parents {
    current = current
      .entry(key.clone())
      .or_insert_with(|| Value::Object(Map::new()))
      .as_object_mut()
      .ok_or_else(|| make_report!("setting `{}` is inside a value that is not a mapping", path.join(".")))?;
  }
  current.insert(last.clone(), value);
  Ok(())
}

fn list_tokens(arg: &Arg, flag: &str, items: &[Value]) -> Option<Vec<String>> {
  if items.is_empty() {
    return None;
  }
  let values = items.iter().map(|item| cli_value(arg, item)).collect_vec();
  if let Some(delimiter) = arg.get_value_delimiter() {
    return Some(vec![flag.to_owned(), values.join(&delimiter.to_string())]);
  }
  let (min, max) = arg
    .get_num_args()
    .map_or((1, 1), |range| (range.min_values(), range.max_values()));
  if max == usize::MAX || (min..=max).contains(&values.len()) {
    return Some([vec![flag.to_owned()], values].concat());
  }
  if max == 0 || values.len() % max != 0 {
    return None;
  }
  Some(
    values
      .chunks(max)
      .flat_map(|chunk| [vec![flag.to_owned()], chunk.to_vec()].concat())
      .collect(),
  )
}

fn cli_value(arg: &Arg, value: &Value) -> String {
  let text = match value {
    Value::String(text) => text.clone(),
    Value::Array(_) | Value::Object(_) => value.to_string(),
    Value::Number(number) => number.to_string(),
    Value::Bool(flag) => flag.to_string(),
    Value::Null => "null".to_owned(),
  };
  arg
    .get_possible_values()
    .into_iter()
    .find(|possible| spelling_key(possible.get_name()) == spelling_key(&text))
    .map_or(text, |possible| possible.get_name().to_owned())
}

fn spelling_key(value: &str) -> String {
  value
    .chars()
    .filter(|c| *c != '-' && *c != '_')
    .flat_map(char::to_lowercase)
    .collect()
}

fn quote_tokens(tokens: &[String]) -> Result<String, Report> {
  tokens
    .iter()
    .map(|token| {
      shlex::try_quote(token)
        .map(|quoted| quoted.into_owned())
        .map_err(|err| make_report!("could not quote `{token}` for the shell: {err}"))
    })
    .collect::<Result<Vec<_>, Report>>()
    .map(|tokens| tokens.join(" "))
}

fn yaml_node(value: &Value) -> Result<Yaml<'static>, Report> {
  let plain = |text: String| Yaml::Representation(Cow::Owned(text), ScalarStyle::Plain, None);
  Ok(match value {
    Value::Null => plain("null".to_owned()),
    Value::Bool(_) | Value::Number(_) => plain(value.to_string()),
    Value::String(_) => {
      let quoted = json_write_str(value, JsonPretty(false))?;
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
