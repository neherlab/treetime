use crate::command::AppCommand;
use crate::config::catalog::{CommandSettings, SettingRole, SettingSpec, command_settings};
use crate::config::settings::setting_ref;
use crate::config::source::{RawDiagnostic, escape_pointer};
use clap::{Arg, ArgAction, Command};
use eyre::Report;
use itertools::Itertools;
use serde_json::{Map, Value};
use treetime_utils::make_report;

pub fn cli_rule_diagnostics(command: AppCommand, config: &Map<String, Value>) -> Result<Vec<RawDiagnostic>, Report> {
  let settings = command_settings(command)?;
  let cli = command.cli_command();
  let given = given_settings(&settings, &cli, config)?;
  let conflicts = given.iter().flat_map(|setting| {
    setting
      .spec
      .conflicts
      .iter()
      .filter_map(|other| given.iter().find(|candidate| candidate.spec.key == *other))
      .filter(|other| setting.spec.key < other.spec.key)
      .map(|other| conflict_diagnostic(setting.spec, other.spec))
  });
  let counts = given.iter().filter_map(value_count_diagnostic);
  Ok(conflicts.chain(counts).collect())
}

struct GivenSetting<'a> {
  spec: &'a SettingSpec,
  arg: &'a Arg,
  value: &'a Value,
}

fn given_settings<'a>(
  settings: &'a CommandSettings,
  cli: &'a Command,
  config: &'a Map<String, Value>,
) -> Result<Vec<GivenSetting<'a>>, Report> {
  settings
    .settings
    .iter()
    .filter(|spec| spec.role != SettingRole::Output)
    .filter_map(|spec| {
      let value = setting_ref(config, &spec.path)?;
      (*value != spec.default_value).then_some((spec, value))
    })
    .map(|(spec, value)| {
      let arg = cli
        .get_arguments()
        .find(|arg| arg.get_long() == spec.flag.strip_prefix("--"))
        .ok_or_else(|| make_report!("setting `{}` has no command-line flag", spec.key))?;
      Ok(is_given(arg, value).then_some(GivenSetting { spec, arg, value }))
    })
    .filter_map(Result::transpose)
    .collect()
}

fn is_given(arg: &Arg, value: &Value) -> bool {
  if !arg.get_action().takes_values() {
    return *value == Value::Bool(true);
  }
  match value {
    Value::Null => false,
    Value::Array(items) => !items.is_empty(),
    _ => true,
  }
}

fn conflict_diagnostic(first: &SettingSpec, second: &SettingSpec) -> RawDiagnostic {
  RawDiagnostic::builder(
    "config::conflict",
    format!("`{}` cannot be used together with `{}`", first.key, second.key),
  )
  .at(setting_pointer(first))
  .key_span(true)
  .help(format!(
    "the command line rejects `{}` together with `{}`; remove one of the two settings",
    first.flag, second.flag
  ))
  .build()
}

fn value_count_diagnostic(setting: &GivenSetting<'_>) -> Option<RawDiagnostic> {
  let Value::Array(items) = setting.value else {
    return None;
  };
  if setting.arg.get_value_delimiter().is_some() {
    return None;
  }
  let range = setting.arg.get_num_args()?;
  let (min, max) = (range.min_values(), range.max_values());
  let repeats = matches!(setting.arg.get_action(), ArgAction::Append);
  let count = items.len();
  let uses = if repeats { count } else { 1 };
  if (1..=uses).any(|times| times * min <= count && count <= times.saturating_mul(max)) {
    return None;
  }
  let takes = if min == max {
    format!("{min}")
  } else if max == usize::MAX {
    format!("at least {min}")
  } else {
    format!("{min} to {max}")
  };
  let each_use = if repeats { " each time it is given" } else { "" };
  Some(
    RawDiagnostic::builder(
      "config::value-count",
      format!(
        "`{}` has {count} {}, but `{}` takes {takes} values{each_use}",
        setting.spec.key,
        plural(count, "value", "values"),
        setting.spec.flag
      ),
    )
    .at(setting_pointer(setting.spec))
    .build(),
  )
}

const fn plural<'a>(count: usize, one: &'a str, many: &'a str) -> &'a str {
  if count == 1 { one } else { many }
}

fn setting_pointer(spec: &SettingSpec) -> String {
  spec
    .path
    .iter()
    .map(|part| format!("/{}", escape_pointer(part)))
    .join("")
}
