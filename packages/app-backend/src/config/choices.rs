use crate::command::AppCommandExt;
use crate::config::catalog::SettingSpec;
use crate::run_checks::SettingPatch;
use app_commands::command::AppCommand;
use app_commands::json_value::JsonValue;
use deser::Serialize;
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use treetime_utils::make_report;

const CHOICES: &[ChoiceDef] = &[
  ChoiceDef {
    name: ChoiceName::ClockRate,
    keys: &["clock_rate", "clock_std_dev"],
    options: &[
      OptionDef {
        name: ChoiceOptionName::Estimate,
        set: &[],
        settings: &[],
      },
      OptionDef {
        name: ChoiceOptionName::Fixed,
        set: &[("clock_rate", PickedValue::Example)],
        settings: &["clock_rate", "clock_std_dev"],
      },
    ],
  },
  ChoiceDef {
    name: ChoiceName::Root,
    keys: &["reroot", "outgroup", "keep_root"],
    options: &[
      OptionDef {
        name: ChoiceOptionName::Reroot,
        set: &[("reroot", PickedValue::MethodUnlessDefault)],
        settings: &["reroot"],
      },
      OptionDef {
        name: ChoiceOptionName::Outgroup,
        set: &[("outgroup", PickedValue::Example)],
        settings: &["outgroup"],
      },
      OptionDef {
        name: ChoiceOptionName::Keep,
        set: &[("keep_root", PickedValue::True)],
        settings: &[],
      },
    ],
  },
];

pub fn setting_choices(command: AppCommand, settings: &[SettingSpec]) -> Result<Vec<SettingChoice>, Report> {
  CHOICES
    .iter()
    .filter(|choice| choice.applies(settings))
    .map(|choice| choice.build(command, settings))
    .try_collect()
}

pub fn active_choices(command: AppCommand, settings: &[SettingSpec], config: &Map<String, Value>) -> Vec<ActiveChoice> {
  let value = |key: &str| config.get(key).filter(|value| !value.is_null());
  let flag = |key: &str| value(key).and_then(Value::as_bool).unwrap_or(false);
  let listed = |key: &str| {
    value(key)
      .and_then(Value::as_array)
      .is_some_and(|items| !items.is_empty())
  };
  CHOICES
    .iter()
    .filter(|choice| choice.applies(settings))
    .map(|choice| ActiveChoice {
      choice: choice.name,
      option: match choice.name {
        ChoiceName::ClockRate => {
          if value("clock_rate").is_some() {
            ChoiceOptionName::Fixed
          } else {
            ChoiceOptionName::Estimate
          }
        },
        ChoiceName::Root => {
          if flag("keep_root") {
            ChoiceOptionName::Keep
          } else if listed("outgroup") {
            ChoiceOptionName::Outgroup
          } else if value("reroot").is_some() || command.reroots_by_default() {
            ChoiceOptionName::Reroot
          } else {
            ChoiceOptionName::Keep
          }
        },
      },
    })
    .collect()
}

/// The option of a choice that a configuration selects.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct ActiveChoice {
  pub choice: ChoiceName,
  pub option: ChoiceOptionName,
}

struct ChoiceDef {
  name: ChoiceName,
  keys: &'static [&'static str],
  options: &'static [OptionDef],
}

impl ChoiceDef {
  fn applies(&self, settings: &[SettingSpec]) -> bool {
    self.keys.iter().all(|key| settings.iter().any(|spec| spec.key == *key))
  }

  fn build(&self, command: AppCommand, settings: &[SettingSpec]) -> Result<SettingChoice, Report> {
    Ok(SettingChoice {
      choice: self.name,
      keys: self.keys.iter().map(|key| (*key).to_owned()).collect(),
      options: self
        .options
        .iter()
        .map(|option| option.build(command, self.keys, settings))
        .try_collect()?,
    })
  }
}

/// Settings that the form shows as one choice between options, for example the clock rate.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct SettingChoice {
  pub choice: ChoiceName,
  /// Keys of the settings the choice controls.
  pub keys: Vec<String>,
  /// The options, in the order the form shows them.
  pub options: Vec<ChoiceOption>,
}

/// A group of settings that the form shows as one choice.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum ChoiceName {
  /// Estimate the clock rate, or fix it.
  ClockRate,
  /// Reroot the tree by a method, root it on an outgroup, or keep the input root.
  Root,
}

struct OptionDef {
  name: ChoiceOptionName,
  set: &'static [(&'static str, PickedValue)],
  settings: &'static [&'static str],
}

impl OptionDef {
  fn build(&self, command: AppCommand, keys: &[&str], settings: &[SettingSpec]) -> Result<ChoiceOption, Report> {
    let patch = keys
      .iter()
      .map(|key| {
        let value = self
          .set
          .iter()
          .find(|(set_key, _)| set_key == key)
          .map(|(_, value)| value.resolve(command, key, settings))
          .transpose()?
          .flatten();
        Ok(SettingPatch {
          key: (*key).to_owned(),
          value,
        })
      })
      .collect::<Result<Vec<_>, Report>>()?;
    Ok(ChoiceOption {
      option: self.name,
      patch,
      settings: self.settings.iter().map(|key| (*key).to_owned()).collect(),
    })
  }
}

/// One option of a choice.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct ChoiceOption {
  pub option: ChoiceOptionName,
  /// Changes that picking the option makes: settings it writes, and settings it removes so that they take their
  /// defaults.
  pub patch: Vec<SettingPatch>,
  /// Keys of the settings the form shows while the option is picked.
  pub settings: Vec<String>,
}

/// An option of a choice.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum ChoiceOptionName {
  /// Estimate the clock rate from the data.
  Estimate,
  /// A fixed clock rate.
  Fixed,
  /// Find the root that fits the molecular clock best.
  Reroot,
  /// Root on the branch above the most recent common ancestor of named tips.
  Outgroup,
  /// Keep the root of the input tree.
  Keep,
}

enum PickedValue {
  True,
  Example,
  MethodUnlessDefault,
}

impl PickedValue {
  fn resolve(&self, command: AppCommand, key: &str, settings: &[SettingSpec]) -> Result<Option<JsonValue>, Report> {
    let spec = || {
      settings
        .iter()
        .find(|spec| spec.key == key)
        .ok_or_else(|| make_report!("setting `{key}` of a choice option is not a setting of {command}"))
    };
    match self {
      Self::True => Ok(Some(JsonValue(Value::Bool(true)))),
      Self::Example => spec()?
        .examples
        .first()
        .cloned()
        .map(Some)
        .ok_or_else(|| make_report!("setting `{key}` has no example value for its choice option")),
      Self::MethodUnlessDefault if command.reroots_by_default() => Ok(None),
      Self::MethodUnlessDefault => spec()?
        .options
        .first()
        .map(|option| Some(JsonValue(Value::String(option.value.clone()))))
        .ok_or_else(|| make_report!("setting `{key}` has no allowed value for its choice option")),
    }
  }
}
