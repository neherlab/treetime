use crate::check_inputs::{InputFacts, InputKind, InputNeed};
use crate::command::AppCommand;
use crate::config::settings::setting_ref;
use crate::config::source::ConfigProblem;
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::iter;
use treetime_utils::vec_of_owned;

const NAMES_SHOWN: usize = 3;

const MONTH_ROUNDING_FACTOR: usize = 4;

const DATE_COLUMN_SETTING: &str = "date_column";

const RULES: [fn(&CheckContext<'_>) -> Vec<RunCheck>; 7] = [
  missing_inputs,
  unreadable_inputs,
  rejected_config,
  tips_without_sequence,
  metadata,
  confidence_without_rate_uncertainty,
  dates_rounded_to_the_month,
];

pub fn run_checks(context: &CheckContext<'_>) -> Vec<RunCheck> {
  RULES.iter().flat_map(|rule| rule(context)).collect()
}

pub fn rejection_messages(rejection: &ConfigRejection<'_>, problems_only: bool) -> Vec<String> {
  let problems = rejection
    .problems
    .iter()
    .map(|problem| match problem.help.as_deref() {
      Some(help) if !help.is_empty() => format!("{} ({help})", problem.message),
      _ => problem.message.clone(),
    })
    .collect_vec();
  if !problems.is_empty() || problems_only {
    return problems;
  }
  vec![
    iter::once(rejection.message)
      .chain(rejection.causes.iter().map(String::as_str))
      .join(": "),
  ]
}

pub struct CheckContext<'a> {
  pub command: AppCommand,
  pub config: Option<&'a Map<String, Value>>,
  pub rejection: Option<ConfigRejection<'a>>,
  pub facts: Option<&'a InputFacts>,
}

pub struct ConfigRejection<'a> {
  pub message: &'a str,
  pub causes: &'a [String],
  pub problems: &'a [ConfigProblem],
}

/// A finding about a configuration and its input files, before a run.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct RunCheck {
  /// Identifier of the check, stable across calls for the same finding.
  pub id: String,
  /// How the finding affects the run.
  pub level: CheckLevel,
  /// The finding, as a sentence.
  pub text: String,
  /// Keys of the settings the finding concerns.
  pub settings: Vec<String>,
  /// Change of settings that resolves the finding, when there is one.
  pub fix: Option<CheckFix>,
}

/// How a finding affects the run.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum CheckLevel {
  /// The run cannot start.
  Block,
  /// The run starts, with a result the user may not expect.
  Warn,
  /// A hint about the inputs.
  Advice,
}

/// Change of settings that resolves a finding.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct CheckFix {
  /// Label of the action, for example `Use covariation`.
  pub label: String,
  /// Settings to set.
  pub patch: Vec<SettingPatch>,
}

/// New value of one setting.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct SettingPatch {
  /// Key path of the setting.
  pub path: Vec<String>,
  /// The value to set.
  pub value: Value,
}

fn missing_inputs(context: &CheckContext<'_>) -> Vec<RunCheck> {
  let Some(config) = context.config else {
    return vec![];
  };
  context
    .command
    .inputs()
    .iter()
    .filter(|input| input.need == InputNeed::Required && !has_path(config.get(input.kind.setting())))
    .map(|input| {
      block(
        format!("missing-{}", input.kind),
        format!("Add {} {} file.", article(input.kind), input.kind),
      )
    })
    .collect()
}

fn unreadable_inputs(context: &CheckContext<'_>) -> Vec<RunCheck> {
  context
    .facts
    .map(|facts| {
      facts
        .problems
        .iter()
        .map(|problem| {
          block(
            format!("unreadable-{}", problem.input),
            format!("The {} cannot be read: {}", problem.input, problem.message),
          )
          .concerning(&[problem.input.setting()])
        })
        .collect()
    })
    .unwrap_or_default()
}

fn rejected_config(context: &CheckContext<'_>) -> Vec<RunCheck> {
  let Some(rejection) = &context.rejection else {
    return vec![];
  };
  let inputs_missing = !missing_inputs(context).is_empty();
  rejection_messages(rejection, inputs_missing)
    .into_iter()
    .enumerate()
    .map(|(index, message)| block(format!("config-{index}"), message))
    .collect()
}

fn tips_without_sequence(context: &CheckContext<'_>) -> Vec<RunCheck> {
  let Some(facts) = context.facts else {
    return vec![];
  };
  let missing = facts.tips_without_sequence.as_deref().unwrap_or_default();
  match &facts.tree {
    Some(tree) if reads(context.command, InputKind::Alignment) && !missing.is_empty() => vec![block(
      "tips-without-sequence",
      format!(
        "{} of {} tree tips have no sequence in the alignment: {}",
        missing.len(),
        tree.tips,
        name_list(missing)
      ),
    )],
    _ => vec![],
  }
}

fn metadata(context: &CheckContext<'_>) -> Vec<RunCheck> {
  let Some(facts) = context.facts else {
    return vec![];
  };
  if !reads(context.command, InputKind::Metadata) {
    return vec![];
  }
  let uses_dates = context.command.uses_dates();
  let mut checks = vec![];

  let missing = facts.tips_without_metadata.as_deref().unwrap_or_default();
  if let Some(tree) = &facts.tree
    && !missing.is_empty()
  {
    let consequence = if uses_dates {
      "get no date"
    } else {
      "get no trait value"
    };
    checks.push(warn(
      "tips-without-metadata",
      format!(
        "{} of {} tree tips have no metadata row and {consequence}: {}",
        missing.len(),
        tree.tips,
        name_list(missing)
      ),
    ));
  }

  if let Some(metadata) = facts.metadata.as_ref().filter(|_| uses_dates) {
    if metadata.date_column.is_none() {
      checks.push(
        block(
          "no-date-column",
          format!(
            "The metadata has no date column. Set the date column to one of: {}.",
            metadata.columns.join(", ")
          ),
        )
        .concerning(&[DATE_COLUMN_SETTING]),
      );
    }
    let unreadable = metadata
      .dates
      .as_ref()
      .map(|dates| dates.unreadable.as_slice())
      .unwrap_or_default();
    if !unreadable.is_empty() {
      checks.push(warn(
        "unreadable-dates",
        format!(
          "{} dates cannot be read and those samples get no date: {}. Use 2015-06-21, 2015-06-XX or 2015.47.",
          unreadable.len(),
          name_list(unreadable)
        ),
      ));
    }
  }
  checks
}

fn confidence_without_rate_uncertainty(context: &CheckContext<'_>) -> Vec<RunCheck> {
  let Some(config) = context.config else {
    return vec![];
  };
  let setting = |key: &str| setting_ref(config, &[key.to_owned()]);
  let confidence = setting("confidence") == Some(&Value::Bool(true));
  let covariation = setting("covariation") == Some(&Value::Bool(true));
  let clock_std_dev = setting("clock_std_dev").is_some_and(|value| !value.is_null());
  if context.command != AppCommand::Timetree || !confidence || covariation || clock_std_dev {
    return vec![];
  }
  vec![RunCheck {
    id: "confidence-without-rate-uncertainty".to_owned(),
    level: CheckLevel::Warn,
    text: "Date intervals need rate uncertainty: without the covariation-aware regression or a clock rate std. dev., this run writes no intervals.".to_owned(),
    settings: vec_of_owned!["confidence", "covariation", "clock_std_dev"],
    fix: Some(CheckFix {
      label: "Use covariation".to_owned(),
      patch: vec![SettingPatch {
        path: vec!["covariation".to_owned()],
        value: Value::Bool(true),
      }],
    }),
  }]
}

fn dates_rounded_to_the_month(context: &CheckContext<'_>) -> Vec<RunCheck> {
  let dates = context
    .facts
    .and_then(|facts| facts.metadata.as_ref())
    .and_then(|metadata| metadata.dates.as_ref());
  match dates {
    Some(dates)
      if context.command.uses_dates()
        && dates.exact_days > 0
        && dates.on_day_1_or_15 * MONTH_ROUNDING_FACTOR > dates.exact_days =>
    {
      vec![RunCheck {
        id: "dates-on-day-1-or-15".to_owned(),
        level: CheckLevel::Advice,
        text: format!(
          "{} of {} dates fall on the 1st or 15th of a month. If only the month is known, write it as 2015-06-XX so TreeTime treats it as a range.",
          dates.on_day_1_or_15, dates.exact_days
        ),
        settings: vec![],
        fix: None,
      }]
    },
    _ => vec![],
  }
}

fn reads(command: AppCommand, kind: InputKind) -> bool {
  command.inputs().iter().any(|input| input.kind == kind)
}

fn has_path(value: Option<&Value>) -> bool {
  match value {
    Some(Value::String(path)) => !path.is_empty(),
    Some(Value::Array(paths)) => paths.iter().any(|path| has_path(Some(path))),
    _ => false,
  }
}

const fn article(kind: InputKind) -> &'static str {
  match kind {
    InputKind::Alignment => "an",
    InputKind::Tree | InputKind::Metadata => "a",
  }
}

fn name_list(names: &[String]) -> String {
  let shown = names.iter().take(NAMES_SHOWN).join(", ");
  if names.len() > NAMES_SHOWN {
    format!("{shown}, and {} more", names.len() - NAMES_SHOWN)
  } else {
    shown
  }
}

fn block(id: impl Into<String>, text: String) -> RunCheck {
  RunCheck {
    id: id.into(),
    level: CheckLevel::Block,
    text,
    settings: vec![],
    fix: None,
  }
}

fn warn(id: impl Into<String>, text: String) -> RunCheck {
  RunCheck {
    id: id.into(),
    level: CheckLevel::Warn,
    text,
    settings: vec![],
    fix: None,
  }
}

impl RunCheck {
  fn concerning(self, settings: &[&str]) -> Self {
    Self {
      settings: settings.iter().map(|setting| (*setting).to_owned()).collect(),
      ..self
    }
  }
}
