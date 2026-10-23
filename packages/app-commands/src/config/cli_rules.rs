use crate::command::AppCommand;
use crate::config::cli_flags::setting_arg;
use crate::config::properties::{PathRole, SettingProperty, setting_properties};
use crate::config::source::{RawDiagnostic, escape_pointer};
use clap::{Arg, ArgAction, Command};
use eyre::Report;
use itertools::{Itertools, izip};
use serde_json::{Map, Value};
use treetime_utils::make_report;

pub fn cli_rule_diagnostics(command: AppCommand, config: &Map<String, Value>) -> Result<Vec<RawDiagnostic>, Report> {
  let cli = command.cli_command();
  let properties = setting_properties(command.config_schema().as_value())?;
  let args: Vec<&Arg> = properties
    .iter()
    .map(|property| setting_arg(&cli, property))
    .try_collect()?;
  let mut given: Vec<GivenSetting<'_>> = izip!(&properties, &args)
    .filter_map(|(property, arg)| {
      let value = config.get(&property.key)?;
      is_given(property, arg, value).then_some((property, *arg, value))
    })
    .map(|(property, arg, value)| GivenSetting::new(property, arg, value))
    .try_collect()?;
  given.sort_by_key(|setting| setting.key);
  let conflicts = given
    .iter()
    .tuple_combinations()
    .filter(|(first, second)| args_conflict(&cli, first.arg, second.arg))
    .map(|(first, second)| conflict_diagnostic(first, second));
  let counts = given.iter().filter_map(value_count_diagnostic);
  Ok(conflicts.chain(counts).collect())
}

pub fn is_given(property: &SettingProperty, arg: &Arg, value: &Value) -> bool {
  if property.path_role == Some(PathRole::Output) {
    return false;
  }
  if !arg.get_action().takes_values() {
    return *value == Value::Bool(true);
  }
  match value {
    Value::Null => false,
    Value::Array(items) => !items.is_empty(),
    _ => true,
  }
}

pub fn setting_conflicts(cli: &Command, properties: &[SettingProperty], args: &[&Arg], arg: &Arg) -> Vec<String> {
  izip!(properties, args)
    .filter(|(_, other)| args_conflict(cli, arg, other))
    .map(|(property, _)| property.key.clone())
    .sorted()
    .dedup()
    .collect()
}

fn args_conflict(cli: &Command, first: &Arg, second: &Arg) -> bool {
  let declares = |arg: &Arg, other: &Arg| {
    cli
      .get_arg_conflicts_with(arg)
      .iter()
      .any(|conflict| conflict.get_id() == other.get_id())
  };
  declares(first, second) || declares(second, first)
}

fn conflict_diagnostic(first: &GivenSetting<'_>, second: &GivenSetting<'_>) -> RawDiagnostic {
  RawDiagnostic::builder(
    "config::conflict",
    format!("`{}` cannot be used together with `{}`", first.key, second.key),
  )
  .at(first.pointer())
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
        setting.key,
        plural(count, "value", "values"),
        setting.flag
      ),
    )
    .at(setting.pointer())
    .build(),
  )
}

struct GivenSetting<'a> {
  key: &'a str,
  flag: String,
  arg: &'a Arg,
  value: &'a Value,
}

impl<'a> GivenSetting<'a> {
  fn new(property: &'a SettingProperty, arg: &'a Arg, value: &'a Value) -> Result<Self, Report> {
    let long = arg
      .get_long()
      .ok_or_else(|| make_report!("config key `{}` has no command-line flag", property.key))?;
    Ok(Self {
      key: &property.key,
      flag: format!("--{long}"),
      arg,
      value,
    })
  }

  fn pointer(&self) -> String {
    format!("/{}", escape_pointer(self.key))
  }
}

const fn plural<'a>(count: usize, one: &'a str, many: &'a str) -> &'a str {
  if count == 1 { one } else { many }
}
