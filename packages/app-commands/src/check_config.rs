use crate::check_inputs::InputFacts;
use crate::command::AppCommand;
use crate::config::catalog::{SettingRole, SettingSpec, command_settings};
use crate::config::choices::{ActiveChoice, active_choices};
use crate::config::code::{ConfigCode, config_code};
use crate::config::resolve_paths::resolve_config_paths_where;
use crate::config::settings::has_path;
use crate::config::source::{ConfigProblem, ConfigSource, InvalidConfig, parse_config_document};
use crate::json_value::SparseConfig;
use crate::run_checks::{CheckContext, ConfigRejection, RunCheck, rejection_messages, run_checks};
use crate::yaml::yaml_text;
use app_datasets::text_schema_command;
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use serde_with::skip_serializing_none;
use std::collections::BTreeSet;
use std::mem;
use std::path::{Path, PathBuf};
use strum::IntoEnumIterator;
use treetime_utils::error::ReportChain;

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
    let text_keys = top_level_keys(&document);
    let (text, document) = with_inputs(&settings, &request.text, document, &request.inputs)?;
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

/// Outcome of checking a configuration without running it.
#[skip_serializing_none]
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "status", rename_all = "kebab-case")]
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

/// Request to check a configuration.
#[skip_serializing_none]
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CheckConfigRequest {
  /// Command the configuration is for, unless the text starts with the `yaml-language-server` schema directive of
  /// another command.
  pub command: AppCommand,
  /// Configuration as YAML or JSON text.
  pub text: String,
  /// Input settings to add when the text does not set them, for example the inputs of a draft that the text is loaded
  /// into.
  #[serde(default)]
  pub inputs: SparseConfig,
  /// Facts about the input files, from `check-inputs`, for the checks that depend on them.
  #[serde(default)]
  pub input_facts: Option<InputFacts>,
  /// Absolute folder that relative paths in the text resolve from: the folder of the config file. Paths in `inputs`
  /// stay as they are. Unset: relative paths resolve from the working directory of the back end.
  #[serde(default)]
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
  let rejection = ConfigRejection {
    message: &message,
    causes: &causes,
    problems: &problems,
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
  resolve_config_paths_where(&mut value, command.config_schema().as_value(), folder, |key_path| {
    key_path.first().is_some_and(|key| text_keys.contains(key))
  })?;
  if let Value::Object(settings) = value {
    *config = settings;
  }
  Ok(())
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
