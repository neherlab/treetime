use crate::command::AppCommandExt;
use crate::config::catalog::{CommandSettings, SettingRole, SettingSpec, command_settings};
use crate::config::effective::value_text;
use crate::results::citation::{Citation, citation};
use crate::run_config::RUN_CONFIG_OUTPUT_DIR;
use app_commands::command::AppCommand;
use app_commands::config::settings::has_path;
use app_commands::yaml::yaml_text;
use app_datasets::schema_directive;
use clap::{Arg, Command};
use deser::Serialize;
use deser::adapters::As;
use deser_serde::Serde;
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use treetime_schema::version_info;
use treetime_utils::io::json::{JsonPretty, json_write_str};
use treetime_utils::make_report;

const OUTPUT_ALL_KEY: &str = "output_all";

const OUTPUT_ALL_FLAG: &str = "--output-all";

const CONFIG_FILE: &str = "run.yaml";

pub fn config_code(command: AppCommand, config: &Map<String, Value>) -> Result<ConfigCode, Report> {
  let settings = command_settings(command)?;
  let cli = command.cli_command();
  let inputs = settings
    .settings
    .iter()
    .filter(|spec| matches!(spec.role, SettingRole::Input | SettingRole::InputTemplate))
    .filter_map(|spec| Some((spec, config.get(&spec.key)?.clone())))
    .filter(|(_, value)| has_path(value))
    .collect_vec();
  let changed = changed_settings(&settings, config);
  let output_dir = config
    .get(OUTPUT_ALL_KEY)
    .and_then(Value::as_str)
    .unwrap_or(RUN_CONFIG_OUTPUT_DIR);

  let command_line = command_line(command, &cli, &inputs, &changed, output_dir)?;
  let yaml = yaml(command, &inputs, &changed, output_dir)?;
  let citation = citation();
  Ok(ConfigCode {
    methods: methods(command, &changed, &citation),
    citation,
    config_file: CONFIG_FILE.to_owned(),
    config_command: config_command(command),
    command_line_text: command_line_text(&command_line),
    yaml_text: format!("{}\n", yaml.iter().map(|line| &line.text).join("\n")),
    command_line,
    yaml,
  })
}

/// A command line and a YAML config that reproduce a configuration.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct ConfigCode {
  /// The command line, one flag per line.
  pub command_line: Vec<CodeLine>,
  /// The command line as a shell reads it, with a line continuation after every flag but the last.
  pub command_line_text: String,
  /// The YAML config, line by line.
  pub yaml: Vec<CodeLine>,
  /// The YAML config as a file holds it.
  pub yaml_text: String,
  /// Name to save the YAML config under.
  pub config_file: String,
  /// The command line that runs the saved YAML config.
  pub config_command: String,
  /// A methods paragraph for a publication: the analysis, the TreeTime version, and the key settings that differ from
  /// their defaults.
  pub methods: String,
  /// The publication to cite for TreeTime.
  pub citation: Citation,
}

fn methods(command: AppCommand, changed: &[(&SettingSpec, Value)], citation: &Citation) -> String {
  let (key, other): (Vec<_>, Vec<_>) = changed
    .iter()
    .partition(|(spec, _)| command.main_settings().contains(&spec.key.as_str()));
  let analysis = format!(
    "We used the {} analysis of TreeTime {} ({}; command `treetime {command}`).",
    lowercase_first(command.title()),
    version_info().version,
    citation.short,
  );
  let key_settings = (!key.is_empty()).then(|| {
    format!(
      "Key settings changed from their defaults: {}.",
      key
        .iter()
        .map(|(spec, value)| format!("{} = {}", lowercase_first(&spec.label), prose_value(value)))
        .join(", ")
    )
  });
  let rest = match (key.is_empty(), other.len()) {
    (true, 0) => "All settings had their default values.".to_owned(),
    (false, 0) => "All other settings had their default values.".to_owned(),
    (_, 1) => "One further setting differs from its default, as the configuration lists it.".to_owned(),
    (_, count) => format!("{count} further settings differ from their defaults, as the configuration lists them."),
  };
  [Some(analysis), key_settings, Some(rest)]
    .into_iter()
    .flatten()
    .join(" ")
}

fn prose_value(value: &Value) -> String {
  match value {
    Value::Bool(true) => "on".to_owned(),
    Value::Bool(false) => "off".to_owned(),
    Value::Null | Value::Number(_) | Value::String(_) | Value::Array(_) | Value::Object(_) => value_text(Some(value)),
  }
}

fn lowercase_first(text: &str) -> String {
  let mut chars = text.chars();
  match (chars.next(), chars.clone().next()) {
    (Some(first), Some(second)) if !second.is_uppercase() => first.to_lowercase().chain(chars).collect(),
    _ => text.to_owned(),
  }
}

fn changed_settings<'a>(settings: &'a CommandSettings, config: &Map<String, Value>) -> Vec<(&'a SettingSpec, Value)> {
  settings
    .settings
    .iter()
    .filter(|spec| spec.role == SettingRole::Setting)
    .filter_map(|spec| {
      let value = config.get(&spec.key)?;
      (Some(value) != spec.default_value.as_deref()).then(|| (spec, value.clone()))
    })
    .collect()
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
            json_write_str(&As::<_, Serde>::new(value), JsonPretty(false))?
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

pub fn setting_tokens(arg: &Arg, flag: &str, value: &Value) -> Option<Vec<String>> {
  if !arg.get_action().takes_values() {
    return (value == &Value::Bool(true)).then(|| vec![flag.to_owned()]);
  }
  match value {
    Value::Null => None,
    Value::Array(items) => list_tokens(arg, flag, items),
    _ => Some(vec![flag.to_owned(), cli_value(value)]),
  }
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
      text: format!("# {}", config_command(command)),
      kind: CodeLineKind::Comment,
    },
  ];
  for (spec, value) in inputs {
    lines.extend(yaml_entry(&spec.key, value, CodeLineKind::Input)?);
  }
  for (spec, value) in changed {
    lines.extend(yaml_entry(&spec.key, value, CodeLineKind::Changed)?);
  }
  lines.extend(yaml_entry(
    OUTPUT_ALL_KEY,
    &Value::String(output_dir.to_owned()),
    CodeLineKind::Output,
  )?);
  Ok(lines)
}

fn config_command(command: AppCommand) -> String {
  format!("treetime {command} --config {CONFIG_FILE}")
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

/// One line of a command line or a YAML config.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct CodeLine {
  /// Text of the line.
  pub text: String,
  /// What the line sets.
  pub kind: CodeLineKind,
}

/// What a line of a command line or a YAML config sets.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
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

fn list_tokens(arg: &Arg, flag: &str, items: &[Value]) -> Option<Vec<String>> {
  if items.is_empty() {
    return None;
  }
  let values = items.iter().map(cli_value).collect_vec();
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

fn cli_value(value: &Value) -> String {
  match value {
    Value::String(text) => text.clone(),
    Value::Array(_) | Value::Object(_) | Value::Number(_) | Value::Bool(_) | Value::Null => value.to_string(),
  }
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
