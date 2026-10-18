use crate::check_inputs::{AlignmentFacts, InputFacts, InputKind, InputNeed};
use crate::command::AppCommandExt;
use app_commands::command::AppCommand;
use app_commands::commands::shared::input_warnings::name_list;
use app_commands::config::settings::{has_path, setting_ref};
use app_commands::config::source::ConfigProblem;
use app_commands::json_value::JsonValue;
use deser::Serialize;
use itertools::Itertools;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use std::iter;
use treetime_schema::skip_serializing_optionals;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::vec_of_owned;

const MONTH_ROUNDING_FACTOR: usize = 4;

const DATE_COLUMN_SETTING: &str = "date_column";

const IGNORE_MISSING_SEQUENCES_SETTING: &str = "ignore_missing_alns";

const RULES: [fn(&CheckContext<'_>) -> Vec<RunCheck>; 9] = [
  missing_inputs,
  unreadable_inputs,
  rejected_config,
  duplicate_names,
  unaligned_sequences,
  tips_without_sequence,
  metadata,
  confidence_without_rate_uncertainty,
  dates_rounded_to_the_month,
];

pub fn run_checks(context: &CheckContext<'_>) -> Vec<RunCheck> {
  RULES.iter().flat_map(|rule| rule(context)).collect()
}

fn unreadable_inputs(context: &CheckContext<'_>) -> Vec<RunCheck> {
  context
    .facts
    .map(|facts| {
      facts
        .problems
        .iter()
        .map(|problem| {
          let id = format!("unreadable-{}", problem.input);
          let id = match problem.partition {
            Some(partition) if several_alignments(context) => format!("{id}-{partition}"),
            _ => id,
          };
          block(id, format!("The {} cannot be read: {}", problem.input, problem.message))
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

fn missing_inputs(context: &CheckContext<'_>) -> Vec<RunCheck> {
  let Some(config) = context.config else {
    return vec![];
  };
  context
    .command
    .inputs()
    .iter()
    .filter(|input| input.need == InputNeed::Required && !config.get(input.kind.setting()).is_some_and(has_path))
    .map(|input| {
      block(
        format!("missing-{}", input.kind),
        format!("Add {} {} file.", article(input.kind), input.kind),
      )
    })
    .collect()
}

fn duplicate_names(context: &CheckContext<'_>) -> Vec<RunCheck> {
  let Some(facts) = context.facts else {
    return vec![];
  };
  let tree = facts
    .tree
    .as_ref()
    .filter(|tree| !tree.duplicate_node_names.is_empty())
    .map(|tree| {
      warn(
        "duplicate-node-names",
        format!(
          "The tree gives the same name to more than one node: {}. Nodes with the same name receive the same data from the other inputs.",
          name_list(&tree.duplicate_node_names)
        ),
      )
      .concerning(&[InputKind::Tree.setting()])
    });
  let alignments = sequence_files(context)
    .filter(|alignment| !alignment.duplicate_names.is_empty())
    .map(|alignment| {
      warn(
        alignment_id(context, "duplicate-sequence-names", alignment),
        format!(
          "The {} has more than one sequence named {}. TreeTime uses the first sequence of each name.",
          alignment_noun(context, alignment),
          name_list(&alignment.duplicate_names)
        ),
      )
      .concerning(&[InputKind::Alignment.setting()])
    });
  tree.into_iter().chain(alignments).collect()
}

fn unaligned_sequences(context: &CheckContext<'_>) -> Vec<RunCheck> {
  sequence_files(context)
    .filter(|alignment| alignment.min_length != alignment.max_length)
    .map(|alignment| {
      block(
        alignment_id(context, "unaligned-sequences", alignment),
        format!(
          "The {} is not aligned: its sequences have lengths {} to {}.",
          alignment_noun(context, alignment),
          alignment.min_length,
          alignment.max_length
        ),
      )
      .concerning(&[InputKind::Alignment.setting()])
    })
    .collect()
}

fn tips_without_sequence(context: &CheckContext<'_>) -> Vec<RunCheck> {
  let Some(tree) = context.facts.and_then(|facts| facts.tree.as_ref()) else {
    return vec![];
  };
  let limit = context.command.missing_sequences();
  let ignored = context
    .config
    .and_then(|config| setting_ref(config, &[IGNORE_MISSING_SEQUENCES_SETTING.to_owned()]))
    == Some(&Value::Bool(true));
  sequence_files(context)
    .filter_map(|alignment| {
      let missing = alignment.tips_without_sequence.as_deref().unwrap_or_default();
      let n_missing = missing.len() + tree.unnamed_tips;
      if n_missing == 0 {
        return None;
      }
      let id = alignment_id(context, "tips-without-sequence", alignment);
      let noun = alignment_noun(context, alignment);
      let over_limit = n_missing.to_f64_exact() > tree.tips.to_f64_exact() / 3.0;
      let text = format!(
        "{n_missing} of {} tree tips have no sequence in the {noun}: {}",
        tree.tips,
        name_list(missing)
      );
      let check = match limit {
        MissingSequences::Limited if over_limit => block(
          id,
          format!("{text}. More than a third of the tips are missing, so the alignment may not belong to the tree."),
        ),
        MissingSequences::LimitedUnlessIgnored if over_limit && !ignored => RunCheck {
          fix: Some(CheckFix {
            label: "Treat missing sequences as unknown".to_owned(),
            patch: vec![SettingPatch {
              path: vec![IGNORE_MISSING_SEQUENCES_SETTING.to_owned()],
              value: Some(JsonValue(Value::Bool(true))),
            }],
          }),
          ..block(
            id,
            format!("{text}. More than a third of the tips are missing, so the alignment may not belong to the tree."),
          )
          .concerning(&[InputKind::Alignment.setting(), IGNORE_MISSING_SEQUENCES_SETTING])
        },
        MissingSequences::Limited | MissingSequences::LimitedUnlessIgnored => {
          warn(id, format!("{text}. TreeTime treats their sequences as unknown."))
        },
        MissingSequences::Unread => return None,
      };
      Some(check)
    })
    .collect()
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
  let Some(config) = context.config.filter(|_| context.command == AppCommand::Timetree) else {
    return vec![];
  };
  let setting = |key: &str| setting_ref(config, &[key.to_owned()]);
  let confidence = setting("confidence") == Some(&Value::Bool(true));
  let covariation = setting("covariation") == Some(&Value::Bool(true));
  let clock_std_dev = setting("clock_std_dev").is_some_and(|value| !value.is_null());
  if !confidence || covariation || clock_std_dev {
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
        value: Some(JsonValue(Value::Bool(true))),
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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MissingSequences {
  Unread,
  Limited,
  LimitedUnlessIgnored,
}

fn sequence_files<'a>(context: &CheckContext<'a>) -> impl Iterator<Item = &'a AlignmentFacts> {
  let read = context.command.missing_sequences() != MissingSequences::Unread;
  context
    .facts
    .filter(|_| read)
    .map(|facts| facts.alignments.as_slice())
    .unwrap_or_default()
    .iter()
}

fn several_alignments(context: &CheckContext<'_>) -> bool {
  context
    .config
    .and_then(|config| config.get(InputKind::Alignment.setting()))
    .and_then(Value::as_array)
    .is_some_and(|paths| paths.len() > 1)
}

fn alignment_id(context: &CheckContext<'_>, id: &str, alignment: &AlignmentFacts) -> String {
  if several_alignments(context) {
    format!("{id}-{}", alignment.partition)
  } else {
    id.to_owned()
  }
}

fn alignment_noun(context: &CheckContext<'_>, alignment: &AlignmentFacts) -> String {
  if several_alignments(context) {
    format!("alignment {} '{}'", alignment.partition, alignment.path.display())
  } else {
    "alignment".to_owned()
  }
}

fn reads(command: AppCommand, kind: InputKind) -> bool {
  command.inputs().iter().any(|input| input.kind == kind)
}

const fn article(kind: InputKind) -> &'static str {
  match kind {
    InputKind::Alignment => "an",
    InputKind::Tree | InputKind::Metadata => "a",
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

/// A finding about a configuration and its input files, before a run.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
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

impl RunCheck {
  fn concerning(self, settings: &[&str]) -> Self {
    Self {
      settings: settings.iter().map(|setting| (*setting).to_owned()).collect(),
      ..self
    }
  }
}

/// How a finding affects the run.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum CheckLevel {
  /// The run cannot start.
  Block,
  /// The run starts, with a result the user may not expect.
  Warn,
  /// A hint about the inputs.
  Advice,
}

/// Change of settings that resolves a finding.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct CheckFix {
  /// Label of the action, for example `Use covariation`.
  pub label: String,
  /// Settings to set.
  pub patch: Vec<SettingPatch>,
}

/// New value of one setting.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct SettingPatch {
  /// Key path of the setting.
  pub path: Vec<String>,
  /// The value to set; absent: the setting is removed and takes its default.
  pub value: Option<JsonValue>,
}
