use crate::check_inputs::InputFacts;
use crate::command::AppCommandExt;
use crate::config::catalog::{SettingRole, SettingSpec, command_settings};
use crate::config::choices::{ActiveChoice, active_choices};
use crate::config::code::{ConfigCode, config_code};
use crate::run_checks::{CheckContext, ConfigRejection, RunCheck, rejection_messages, run_checks};
use app_commands::command::AppCommand;
use app_commands::commands::shared::argument_error::ArgumentError;
use app_commands::config::resolve_paths::resolve_config_paths_where;
use app_commands::config::settings::has_path;
use app_commands::config::source::{ConfigProblem, ConfigSource, InvalidConfig, parse_config_document};
use app_commands::json_value::SparseConfig;
use app_commands::yaml::{yaml_document, yaml_text};
use app_datasets::text_schema_command;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use std::collections::BTreeSet;
use std::mem;
use std::path::{Path, PathBuf};
use strum::IntoEnumIterator;
use treetime_schema::field_default;
use treetime_schema::skip_serializing_optionals;
use treetime_utils::error::ReportChain;
use treetime_utils::make_error;

const SOURCE_NAME: &str = "config.yaml";

pub fn check_config(request: &CheckConfigRequest) -> CheckConfigResponse {
  let command = text_schema_command(&request.text)
    .and_then(|name| AppCommand::iter().find(|candidate| <&str>::from(*candidate) == name))
    .unwrap_or(request.command);
  let facts = request.input_facts.as_ref();
  let source = ConfigSource::new(SOURCE_NAME, request.text.as_str());
  let settings = match command_settings(command) {
    Ok(settings) => settings.settings,
    Err(report) => return invalid(command, &report, None, &[], facts),
  };
  let prepared = parse_config_document(&source, &request.text).and_then(|document| {
    let (text, document) = with_replaced_inputs(&settings, &request.text, document, &request.replace_inputs)?;
    let mut text_keys = top_level_keys(&document);
    text_keys.retain(|key| !request.replace_inputs.contains_key(key));
    let (text, document) = with_inputs(&settings, &text, document, &request.inputs)?;
    Ok((text, text_keys, command.config_over_defaults(&document)?))
  });
  let (text, text_keys, config) = match prepared {
    Ok(prepared) => prepared,
    Err(report) => return invalid(command, &report, None, &settings, facts),
  };
  let resolved = command.prepare_text(SOURCE_NAME, &text).and_then(|prepared| {
    let mut config = prepared.config;
    if let Some(folder) = &request.folder {
      resolve_text_paths(command, &mut config, folder, &text_keys)?;
    }
    command.remove_output_paths(&mut config)?;
    Ok((config_code(command, &config)?, config))
  });
  match resolved {
    Ok((code, config)) => CheckConfigResponse::Valid {
      command,
      choices: active_choices(&settings, &config),
      checks: run_checks(&CheckContext {
        command,
        config: Some(&config),
        rejection: None,
        facts,
      }),
      config: SparseConfig(config),
      code,
    },
    Err(report) => invalid(command, &report, Some(&config), &settings, facts),
  }
}

/// Request to check a configuration.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct CheckConfigRequest {
  /// Command the configuration is for, unless the text starts with the `yaml-language-server` schema directive of
  /// another command.
  pub command: AppCommand,
  /// Configuration as YAML or JSON text.
  pub text: String,
  /// Input settings to add when the text does not set them, for example the inputs of a draft that the text is loaded
  /// into.
  #[schemars(default, transform = field_default::<SparseConfig>)]
  #[deser(default)]
  pub inputs: SparseConfig,
  /// Input settings whose values replace the paths that the text sets, for example the uploaded copies of the files
  /// that an example configuration reads. The replaced values pass the same checks as the text, and paths in
  /// `replace_inputs` stay as they are. With replacements, the text is checked as the YAML document with the
  /// replaced values and without its comments, and the places of problems refer to that document.
  #[schemars(default, transform = field_default::<SparseConfig>)]
  #[deser(default)]
  pub replace_inputs: SparseConfig,
  /// Facts about the input files, from `check-inputs`, for the checks that depend on them.
  #[schemars(default)]
  #[deser(default)]
  pub input_facts: Option<InputFacts>,
  /// Absolute folder that relative paths in the text resolve from: the folder of the config file. Paths in `inputs`
  /// stay as they are. Unset: relative paths resolve from the working directory of the back end.
  #[schemars(default)]
  #[deser(default)]
  pub folder: Option<PathBuf>,
}

fn invalid(
  command: AppCommand,
  report: &Report,
  config: Option<&Map<String, Value>>,
  settings: &[SettingSpec],
  facts: Option<&InputFacts>,
) -> CheckConfigResponse {
  let invalid = report.downcast_ref::<InvalidConfig>();
  let ReportChain { message, causes } = ReportChain::of(report);
  let problems = InvalidConfig::problems_of(report);
  let argument_settings = report
    .downcast_ref::<ArgumentError>()
    .map(|error| error.settings.as_slice())
    .unwrap_or_default();
  let rejection = ConfigRejection {
    message: &message,
    causes: &causes,
    problems: &problems,
    settings: argument_settings,
  };
  let messages = rejection_messages(&rejection, false);
  let checks = run_checks(&CheckContext {
    command,
    config,
    rejection: Some(rejection),
    facts,
  });
  CheckConfigResponse::Invalid {
    command,
    choices: config
      .map(|config| active_choices(settings, config))
      .unwrap_or_default(),
    rendered: invalid.map(|invalid| invalid.rendered.clone()),
    messages,
    checks,
    message,
    causes,
    problems,
  }
}

/// Outcome of checking a configuration without running it.
#[derive(Clone, Debug, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(tag = "status", rename_all = "kebab-case")]
#[deser(tag = "status", rename_all = "kebab-case")]
pub enum CheckConfigResponse {
  /// The configuration is accepted.
  Valid {
    /// Command the configuration is for.
    command: AppCommand,
    /// The configuration with every default filled in and without output paths, which the app sets for each run.
    config: SparseConfig,
    /// The command line and the YAML config that reproduce the configuration.
    code: ConfigCode,
    /// Findings about the configuration and its input files.
    checks: Vec<RunCheck>,
    /// The option of each choice of the command that the configuration selects.
    choices: Vec<ActiveChoice>,
  },
  /// The configuration is rejected.
  Invalid {
    /// Command the configuration is for.
    command: AppCommand,
    /// The error, as the CLI prints it.
    message: String,
    /// The errors that caused `message`, outermost first.
    causes: Vec<String>,
    /// Problems found by parsing and by the schema check, empty for other errors.
    problems: Vec<ConfigProblem>,
    /// The problems drawn against the configuration text, when the text could be parsed.
    rendered: Option<String>,
    /// The problems as a user reads them: each parse and schema problem with its help, or the error and its causes.
    messages: Vec<String>,
    /// Findings about the configuration and its input files; the rejection is among them.
    checks: Vec<RunCheck>,
    /// The option of each choice of the command that the configuration selects, when the configuration could be
    /// merged over the defaults.
    choices: Vec<ActiveChoice>,
  },
}

fn top_level_keys(document: &Value) -> BTreeSet<String> {
  document
    .as_object()
    .map(|settings| settings.keys().cloned().collect())
    .unwrap_or_default()
}

fn resolve_text_paths(
  command: AppCommand,
  config: &mut Map<String, Value>,
  folder: &Path,
  text_keys: &BTreeSet<String>,
) -> Result<(), Report> {
  let mut value = Value::Object(mem::take(config));
  resolve_config_paths_where(&mut value, command.config_schema().as_value(), folder, |key| {
    text_keys.contains(key)
  })?;
  if let Value::Object(settings) = value {
    *config = settings;
  }
  Ok(())
}

fn with_replaced_inputs(
  settings: &[SettingSpec],
  text: &str,
  document: Value,
  replacements: &Map<String, Value>,
) -> Result<(String, Value), Report> {
  if replacements.is_empty() {
    return Ok((text.to_owned(), document));
  }
  if let Some(key) = replacements.keys().find(|key| {
    !settings
      .iter()
      .any(|spec| spec.role == SettingRole::Input && spec.key == **key)
  }) {
    return make_error!("`replace_inputs` names `{key}`, which is not an input setting of the command");
  }
  let mut document_settings = match document {
    Value::Object(document_settings) => document_settings,
    Value::Null => Map::new(),
    _ => return make_error!("a command configuration must be a mapping of settings"),
  };
  for (key, value) in replacements {
    document_settings.insert(key.clone(), value.clone());
  }
  let document = Value::Object(document_settings);
  Ok((yaml_document(&document)?, document))
}

fn with_inputs(
  settings: &[SettingSpec],
  text: &str,
  mut document: Value,
  inputs: &Map<String, Value>,
) -> Result<(String, Value), Report> {
  let Value::Object(document_settings) = &mut document else {
    return Ok((text.to_owned(), document));
  };
  let mut added = String::new();
  for spec in settings {
    if spec.role != SettingRole::Input || document_settings.contains_key(&spec.key) {
      continue;
    }
    if let Some(value) = inputs.get(&spec.key).filter(|value| has_path(value)) {
      added.push_str(&yaml_text(&spec.key, value)?);
      added.push('\n');
      document_settings.insert(spec.key.clone(), value.clone());
    }
  }
  if added.is_empty() {
    return Ok((text.to_owned(), document));
  }
  let separator = if text.is_empty() || text.ends_with('\n') {
    ""
  } else {
    "\n"
  };
  Ok((format!("{text}{separator}{added}"), document))
}
