use crate::config::catalog::SettingSpec;
use crate::config::choices::{ChoiceName, SettingChoice};
use app_commands::command::AppCommand;
use app_commands::commands::timetree::coalescent::TimescaleUnit;
use deser::Serialize;
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use treetime_schema::skip_serializing_optionals;
use treetime_utils::{make_error, make_report};

const ROWS: [RowDef; 8] = [
  RowDef {
    label: "Clock rate",
    description: Some("Substitutions per site per year"),
    settings: &[("clock_rate", "Rate"), ("clock_std_dev", "Std. dev.")],
    control: ControlDef::Choice(ChoiceName::ClockRate),
  },
  RowDef {
    label: "Confidence intervals",
    description: Some("Marginal dating uncertainty of nodes"),
    settings: &[("confidence", "Confidence intervals")],
    control: ControlDef::Switches(&[SwitchDef {
      key: "confidence",
      text: "Compute confidence intervals",
      inverted: false,
    }]),
  },
  RowDef {
    label: "Coalescent prior",
    description: Some("Tree prior on node times"),
    settings: &[
      ("coalescent", "Time scale"),
      ("n_skyline", "Skyline grid points"),
      ("skyline_stiffness", "Skyline stiffness"),
    ],
    control: ControlDef::Coalescent {
      setting: "coalescent",
      skyline_settings: &["n_skyline", "skyline_stiffness"],
      options: &[
        OptionDef {
          value: "none",
          label: "None",
          note: Some("The node times have no coalescent prior."),
        },
        OptionDef {
          value: "fixed",
          label: "Fixed constant",
          note: Some("Constant population size with a time scale in years."),
        },
        OptionDef {
          value: "opt",
          label: "Optimized constant",
          note: Some("Constant population size with a time scale fitted to the tree."),
        },
        OptionDef {
          value: "skyline",
          label: "Skyline",
          note: Some("Piecewise-constant population size on a grid of time points."),
        },
      ],
    },
  },
  RowDef {
    label: "Root",
    description: Some("Root position of the input tree"),
    settings: &[("reroot", "Method"), ("outgroup", "Tips"), ("keep_root", "Keep root")],
    control: ControlDef::Choice(ChoiceName::Root),
  },
  RowDef {
    label: "Clock filter",
    description: Some("Outlier threshold in interquartile ranges"),
    settings: &[("clock_filter", "Threshold")],
    control: ControlDef::ClockFilter {
      setting: "clock_filter",
      off: "Off: all samples stay in the regression.",
      on: "Samples farther from the regression than the threshold, in interquartile ranges, are ignored.",
    },
  },
  RowDef {
    label: "Polytomies",
    description: Some("Nodes with more than two children"),
    settings: &[("keep_polytomies", "Keep polytomies")],
    control: ControlDef::Switches(&[SwitchDef {
      key: "keep_polytomies",
      text: "Keep polytomies unresolved",
      inverted: false,
    }]),
  },
  RowDef {
    label: "Substitution model",
    description: None,
    settings: &[("model", "Substitution model"), ("model_params", "Parameters")],
    control: ControlDef::SubstitutionModel {
      model_setting: "model",
      parameters_setting: "model_params",
      parameters_example: "kappa=0.2 pis=0.25,0.25,0.25,0.25",
    },
  },
  RowDef {
    label: "Topology optimizations",
    description: Some("On by default; switch off to keep the topology as it is"),
    settings: &[
      ("no_merge_shared_mutations", "Merge shared mutations"),
      ("no_collapse_short_branches", "Collapse short branches"),
      ("no_hoist_reversions", "Hoist reversions"),
    ],
    control: ControlDef::Switches(&[
      SwitchDef {
        key: "no_merge_shared_mutations",
        text: "Merge polytomy siblings that share mutations",
        inverted: true,
      },
      SwitchDef {
        key: "no_collapse_short_branches",
        text: "Collapse short internal branches",
        inverted: true,
      },
      SwitchDef {
        key: "no_hoist_reversions",
        text: "Hoist reversions in polytomies",
        inverted: true,
      },
    ]),
  },
];

pub fn main_setting_rows(
  command: AppCommand,
  main_settings: &[&str],
  settings: &[SettingSpec],
  choices: &[SettingChoice],
) -> Result<Vec<MainSettingRow>, Report> {
  main_settings
    .iter()
    .map(|key| {
      ROWS
        .iter()
        .position(|row| row.has(key))
        .map_or(RowRef::Setting(key), RowRef::Row)
    })
    .unique()
    .map(|row| match row {
      RowRef::Row(index) => ROWS[index].build(command, main_settings, choices),
      RowRef::Setting(key) => {
        let spec = settings
          .iter()
          .find(|spec| spec.key == key)
          .ok_or_else(|| make_report!("main setting `{key}` is not a setting of {command}"))?;
        Ok(MainSettingRow {
          label: spec.label.clone(),
          description: None,
          settings: vec![RowSetting {
            key: spec.key.clone(),
            label: spec.label.clone(),
          }],
          control: RowControl::Setting,
        })
      },
    })
    .try_collect()
}

/// One row of the main settings: one or more settings that the form shows and edits together.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct MainSettingRow {
  /// Name of the row, for example `Coalescent prior`.
  pub label: String,
  /// Short line under the name, for example `Tree prior on node times`.
  pub description: Option<String>,
  /// The settings of the row, in the order the row shows them.
  pub settings: Vec<RowSetting>,
  /// The control that edits the settings.
  pub control: RowControl,
}

/// One setting of a row of the main settings.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct RowSetting {
  /// Key of the setting.
  pub key: String,
  /// Name of the setting inside the row, for example `Rate` in the row `Clock rate`.
  pub label: String,
}

/// Control of a row of the main settings.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[schemars(tag = "kind", rename_all = "kebab-case")]
#[deser(tag = "kind", rename_all = "kebab-case")]
pub enum RowControl {
  /// A choice between options; the choice of the command with this name holds the options and their texts.
  Choice { choice: ChoiceName },
  /// Switches, one per setting.
  Switches { switches: Vec<RowSwitch> },
  /// The coalescent prior: no prior, a fixed time scale in years, an optimized time scale, or a skyline.
  Coalescent {
    /// Key of the setting that holds the prior: an option name, or a fixed time scale.
    setting: String,
    /// Keys of the settings of the skyline, which apply only while the skyline is picked.
    skyline_settings: Vec<String>,
    /// Unit suffix of a fixed time scale in years, as the setting takes it after the number: `y` in `50y`.
    years_suffix: String,
    options: Vec<RowOption>,
  },
  /// The clock filter threshold.
  ClockFilter {
    /// Key of the setting that holds the threshold.
    setting: String,
    /// Note while the threshold is 0.
    off: String,
    /// Note while the threshold is above 0.
    on: String,
  },
  /// The substitution model and its parameters.
  SubstitutionModel {
    /// Key of the setting that holds the model.
    model_setting: String,
    /// Key of the setting that holds the parameters of the model.
    parameters_setting: String,
    /// Example of the parameters, as the form shows them in the empty field.
    parameters_example: String,
  },
  /// The control of the one setting of the row.
  Setting,
}

/// A switch of a row of the main settings.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct RowSwitch {
  /// Key of the setting the switch turns on or off.
  pub key: String,
  /// Text beside the switch.
  pub text: String,
  /// Whether the switch is on while the setting is off: the setting turns off a step that runs by default.
  pub inverted: bool,
}

/// An option of a row of the main settings.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct RowOption {
  /// The option, for example `skyline`.
  pub value: String,
  /// Name of the option, for example `Skyline`.
  pub label: String,
  /// What the option does, as the form shows it while the option is picked.
  pub note: Option<String>,
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum RowRef<'a> {
  Row(usize),
  Setting(&'a str),
}

struct RowDef {
  label: &'static str,
  description: Option<&'static str>,
  settings: &'static [(&'static str, &'static str)],
  control: ControlDef,
}

impl RowDef {
  fn has(&self, key: &str) -> bool {
    self.settings.iter().any(|(member, _)| *member == key)
  }

  fn build(
    &self,
    command: AppCommand,
    main_settings: &[&str],
    choices: &[SettingChoice],
  ) -> Result<MainSettingRow, Report> {
    let shown = |key: &str| main_settings.contains(&key);
    let control = match self.control {
      ControlDef::Choice(name) => {
        if !choices.iter().any(|choice| choice.choice == name) {
          return make_error!("the row `{}` of {command} has no choice `{name:?}`", self.label);
        }
        RowControl::Choice { choice: name }
      },
      ControlDef::Switches(switches) => RowControl::Switches {
        switches: switches
          .iter()
          .filter(|switch| shown(switch.key))
          .map(|switch| RowSwitch {
            key: switch.key.to_owned(),
            text: switch.text.to_owned(),
            inverted: switch.inverted,
          })
          .collect(),
      },
      ControlDef::Coalescent {
        setting,
        skyline_settings,
        options,
      } => RowControl::Coalescent {
        setting: setting.to_owned(),
        skyline_settings: skyline_settings.iter().map(|key| (*key).to_owned()).collect(),
        years_suffix: TimescaleUnit::Years.to_string(),
        options: options
          .iter()
          .map(|option| RowOption {
            value: option.value.to_owned(),
            label: option.label.to_owned(),
            note: option.note.map(str::to_owned),
          })
          .collect(),
      },
      ControlDef::ClockFilter { setting, off, on } => RowControl::ClockFilter {
        setting: setting.to_owned(),
        off: off.to_owned(),
        on: on.to_owned(),
      },
      ControlDef::SubstitutionModel {
        model_setting,
        parameters_setting,
        parameters_example,
      } => RowControl::SubstitutionModel {
        model_setting: model_setting.to_owned(),
        parameters_setting: parameters_setting.to_owned(),
        parameters_example: parameters_example.to_owned(),
      },
    };
    Ok(MainSettingRow {
      label: self.label.to_owned(),
      description: self.description.map(str::to_owned),
      settings: self
        .settings
        .iter()
        .filter(|(key, _)| shown(key))
        .map(|(key, label)| RowSetting {
          key: (*key).to_owned(),
          label: (*label).to_owned(),
        })
        .collect(),
      control,
    })
  }
}

enum ControlDef {
  Choice(ChoiceName),
  Switches(&'static [SwitchDef]),
  Coalescent {
    setting: &'static str,
    skyline_settings: &'static [&'static str],
    options: &'static [OptionDef],
  },
  ClockFilter {
    setting: &'static str,
    off: &'static str,
    on: &'static str,
  },
  SubstitutionModel {
    model_setting: &'static str,
    parameters_setting: &'static str,
    parameters_example: &'static str,
  },
}

struct SwitchDef {
  key: &'static str,
  text: &'static str,
  inverted: bool,
}

struct OptionDef {
  value: &'static str,
  label: &'static str,
  note: Option<&'static str>,
}
