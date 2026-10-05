use crate::config::catalog::SettingSpec;
use crate::config::settings::setting_ref;
use crate::json_value::JsonValue;
use crate::run_checks::SettingPatch;
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use treetime::timetree::coalescent_timescale::{CoalescentMode, coalescent_mode};
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
    name: ChoiceName::CoalescentPrior,
    keys: &[
      "coalescent",
      "coalescent_opt",
      "coalescent_skyline",
      "skyline_n_points",
      "skyline_stiffness",
    ],
    options: &[
      OptionDef {
        name: ChoiceOptionName::None,
        set: &[],
        settings: &[],
      },
      OptionDef {
        name: ChoiceOptionName::Fixed,
        set: &[("coalescent", PickedValue::Example)],
        settings: &["coalescent"],
      },
      OptionDef {
        name: ChoiceOptionName::Optimized,
        set: &[("coalescent_opt", PickedValue::True)],
        settings: &[],
      },
      OptionDef {
        name: ChoiceOptionName::Skyline,
        set: &[("coalescent_skyline", PickedValue::True)],
        settings: &["skyline_n_points", "skyline_stiffness"],
      },
    ],
  },
  ChoiceDef {
    name: ChoiceName::Root,
    keys: &["reroot", "keep_root"],
    options: &[
      OptionDef {
        name: ChoiceOptionName::Reroot,
        set: &[],
        settings: &["reroot"],
      },
      OptionDef {
        name: ChoiceOptionName::Keep,
        set: &[("keep_root", PickedValue::True)],
        settings: &[],
      },
    ],
  },
];

pub fn setting_choices(settings: &[SettingSpec]) -> Result<Vec<SettingChoice>, Report> {
  CHOICES
    .iter()
    .filter(|choice| choice.applies(settings))
    .map(|choice| choice.build(settings))
    .try_collect()
}

pub fn active_choices(settings: &[SettingSpec], config: &Map<String, Value>) -> Vec<ActiveChoice> {
  let value = |key: &str| setting_ref(config, &[key.to_owned()]);
  let flag = |key: &str| value(key).and_then(Value::as_bool).unwrap_or(false);
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
        ChoiceName::CoalescentPrior => {
          match coalescent_mode(
            value("coalescent").and_then(Value::as_f64),
            flag("coalescent_opt"),
            flag("coalescent_skyline"),
          ) {
            CoalescentMode::Disabled => ChoiceOptionName::None,
            CoalescentMode::Fixed(_) => ChoiceOptionName::Fixed,
            CoalescentMode::Constant => ChoiceOptionName::Optimized,
            CoalescentMode::Skyline => ChoiceOptionName::Skyline,
          }
        },
        ChoiceName::Root => {
          if flag("keep_root") {
            ChoiceOptionName::Keep
          } else {
            ChoiceOptionName::Reroot
          }
        },
      },
    })
    .collect()
}

/// Settings that the form shows as one choice between options, for example the coalescent prior.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct SettingChoice {
  pub choice: ChoiceName,
  /// Keys of the settings the choice controls.
  pub keys: Vec<String>,
  /// The options, in the order the form shows them.
  pub options: Vec<ChoiceOption>,
}

/// One option of a choice.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct ChoiceOption {
  pub option: ChoiceOptionName,
  /// Changes that picking the option makes: settings it writes, and settings it removes so that they take their
  /// defaults.
  pub patch: Vec<SettingPatch>,
  /// Keys of the settings the form shows while the option is picked.
  pub settings: Vec<String>,
}

/// The option of a choice that a configuration selects.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct ActiveChoice {
  pub choice: ChoiceName,
  pub option: ChoiceOptionName,
}

/// A group of settings that the form shows as one choice.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum ChoiceName {
  /// Estimate the clock rate, or fix it.
  ClockRate,
  /// The prior on node times.
  CoalescentPrior,
  /// Reroot the tree, or keep the input root.
  Root,
}

/// An option of a choice.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum ChoiceOptionName {
  /// Estimate the clock rate from the data.
  Estimate,
  /// A fixed value: the clock rate, or the coalescent time scale Tc.
  Fixed,
  /// No coalescent prior.
  None,
  /// One constant coalescent time scale Tc, fitted to the tree.
  Optimized,
  /// A coalescent time scale that changes over time.
  Skyline,
  /// Find the root that fits the molecular clock best.
  Reroot,
  /// Keep the root of the input tree.
  Keep,
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

  fn build(&self, settings: &[SettingSpec]) -> Result<SettingChoice, Report> {
    Ok(SettingChoice {
      choice: self.name,
      keys: self.keys.iter().map(|key| (*key).to_owned()).collect(),
      options: self
        .options
        .iter()
        .map(|option| option.build(self.keys, settings))
        .try_collect()?,
    })
  }
}

struct OptionDef {
  name: ChoiceOptionName,
  set: &'static [(&'static str, PickedValue)],
  settings: &'static [&'static str],
}

impl OptionDef {
  fn build(&self, keys: &[&str], settings: &[SettingSpec]) -> Result<ChoiceOption, Report> {
    let patch = keys
      .iter()
      .map(|key| {
        let value = self
          .set
          .iter()
          .find(|(set_key, _)| set_key == key)
          .map(|(_, value)| value.resolve(key, settings))
          .transpose()?;
        Ok(SettingPatch {
          path: vec![(*key).to_owned()],
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

enum PickedValue {
  True,
  Example,
}

impl PickedValue {
  fn resolve(&self, key: &str, settings: &[SettingSpec]) -> Result<JsonValue, Report> {
    match self {
      Self::True => Ok(JsonValue(Value::Bool(true))),
      Self::Example => settings
        .iter()
        .find(|spec| spec.key == key)
        .and_then(|spec| spec.examples.first().cloned())
        .ok_or_else(|| make_report!("setting `{key}` has no example value for its choice option")),
    }
  }
}
