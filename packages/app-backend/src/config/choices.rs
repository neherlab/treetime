use crate::command::AppCommandExt;
use crate::config::catalog::SettingSpec;
use crate::run_checks::SettingPatch;
use app_commands::command::AppCommand;
use app_commands::json_value::JsonValue;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use treetime_schema::skip_serializing_optionals;
use treetime_utils::make_report;

const CHOICES: &[ChoiceDef] = &[
  ChoiceDef {
    name: ChoiceName::ClockRate,
    keys: &["clock_rate", "clock_std_dev"],
    options: &[
      OptionDef {
        name: ChoiceOptionName::Estimate,
        label: "Estimate from data",
        note: Some("Estimated by root-to-tip regression, then refined during dating."),
        set: &[],
        settings: &[],
      },
      OptionDef {
        name: ChoiceOptionName::Fixed,
        label: "Fixed rate",
        note: None,
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
        label: "Reroot",
        note: None,
        set: &[("reroot", PickedValue::MethodUnlessDefault)],
        settings: &["reroot"],
      },
      OptionDef {
        name: ChoiceOptionName::Outgroup,
        label: "Outgroup",
        note: Some("The root is placed on the branch above the most recent common ancestor of the named tips."),
        set: &[],
        settings: &["outgroup"],
      },
      OptionDef {
        name: ChoiceOptionName::Keep,
        label: "Keep root",
        note: Some("The input root is kept; use this when the input tree is already rooted."),
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

pub fn active_choices(
  command: AppCommand,
  choices: &[SettingChoice],
  settings: &[SettingSpec],
  config: &Map<String, Value>,
  picked: &[ActiveChoice],
) -> Vec<ActiveChoice> {
  choices
    .iter()
    .map(|choice| ActiveChoice {
      choice: choice.choice,
      option: picked
        .iter()
        .find(|active| active.choice == choice.choice)
        .and_then(|active| choice.options.iter().find(|option| option.option == active.option))
        .filter(|option| option.agrees(settings, config))
        .map_or_else(
          || selected_option(choice.choice, command, config),
          |option| option.option,
        ),
    })
    .collect()
}

/// The option of a choice that a configuration selects, or that the form picked.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
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
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
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
  label: &'static str,
  note: Option<&'static str>,
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
      label: self.label.to_owned(),
      note: self.note.map(str::to_owned),
      patch,
      settings: self.settings.iter().map(|key| (*key).to_owned()).collect(),
    })
  }
}

/// One option of a choice.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct ChoiceOption {
  pub option: ChoiceOptionName,
  /// Name of the option, for example `Fixed rate`.
  pub label: String,
  /// What the option does, as the form shows it while the option is picked.
  pub note: Option<String>,
  /// Changes that picking the option makes: settings it writes, and settings it removes so that they take their
  /// defaults.
  pub patch: Vec<SettingPatch>,
  /// Keys of the settings the form shows while the option is picked.
  pub settings: Vec<String>,
}

impl ChoiceOption {
  fn agrees(&self, settings: &[SettingSpec], config: &Map<String, Value>) -> bool {
    self.patch.iter().all(|patch| {
      let expected = patch
        .value
        .as_ref()
        .or_else(|| {
          settings
            .iter()
            .find(|spec| spec.key == patch.key)
            .and_then(|spec| spec.default_value.as_ref())
        })
        .map(|value| &value.0)
        .filter(|value| !value.is_null());
      expected == config.get(&patch.key).filter(|value| !value.is_null())
    })
  }
}

/// An option of a choice.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
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

fn selected_option(choice: ChoiceName, command: AppCommand, config: &Map<String, Value>) -> ChoiceOptionName {
  let value = |key: &str| config.get(key).filter(|value| !value.is_null());
  let flag = |key: &str| value(key).and_then(Value::as_bool).unwrap_or(false);
  let listed = |key: &str| {
    value(key)
      .and_then(Value::as_array)
      .is_some_and(|items| !items.is_empty())
  };
  match choice {
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
  }
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
