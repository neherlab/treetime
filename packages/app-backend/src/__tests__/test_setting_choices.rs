#[cfg(test)]
mod tests {
  use crate::__tests__::test_check_config::tests::helpers::fresh_draft;
  use crate::check_config::CheckConfigResponse;
  use crate::config::catalog::command_settings;
  use crate::config::choices::{ActiveChoice, ChoiceName, ChoiceOption, ChoiceOptionName, SettingChoice};
  use crate::run_checks::SettingPatch;
  use app_commands::command::AppCommand;
  use app_commands::config::properties::setting_properties;
  use app_commands::json_value::JsonValue;
  use helpers::{active_choices_of, check};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use std::collections::BTreeSet;
  use strum::IntoEnumIterator;
  use treetime_utils::{o, vec_of_owned};

  #[test]
  fn test_choices_name_only_settings_of_their_command() {
    let unknown = AppCommand::iter()
      .flat_map(|command| {
        let settings = command_settings(command).unwrap();
        let keys: BTreeSet<String> = setting_properties(command.config_schema().as_value())
          .unwrap()
          .into_iter()
          .map(|setting| setting.key)
          .collect();
        let named = settings
          .choices
          .iter()
          .flat_map(|choice| {
            choice
              .keys
              .iter()
              .cloned()
              .chain(choice.options.iter().flat_map(|option| {
                option
                  .settings
                  .iter()
                  .cloned()
                  .chain(option.patch.iter().map(|patch| patch.key.clone()))
              }))
          })
          .chain(settings.main_settings.iter().cloned())
          .collect_vec();
        named
          .into_iter()
          .filter(|key| !keys.contains(key))
          .map(|key| format!("{command}: {key}"))
          .collect_vec()
      })
      .collect_vec();
    assert_eq!(Vec::<String>::new(), unknown);
  }

  #[test]
  fn test_choices_of_each_command() {
    let choices = AppCommand::iter()
      .map(|command| {
        let names = command_settings(command)
          .unwrap()
          .choices
          .iter()
          .map(|choice| choice.choice)
          .collect_vec();
        (command, names)
      })
      .collect_vec();
    assert_eq!(
      vec![
        (AppCommand::Timetree, vec![ChoiceName::ClockRate, ChoiceName::Root]),
        (AppCommand::Clock, vec![ChoiceName::Root]),
        (AppCommand::Ancestral, vec![]),
        (AppCommand::Homoplasy, vec![]),
        (AppCommand::Mugration, vec![]),
        (AppCommand::Optimize, vec![ChoiceName::Root]),
        (AppCommand::Prune, vec![]),
      ],
      choices
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree(AppCommand::Timetree, None)]
  #[case::clock(   AppCommand::Clock,    None)]
  #[case::optimize(AppCommand::Optimize, Some(json!("min-dev")))]
  #[trace]
  fn test_choices_root_picks_a_method_an_outgroup_or_the_input_root(
    #[case] command: AppCommand,
    #[case] reroot_method: Option<Value>,
  ) {
    let patch = |reroot: Option<Value>, outgroup: Option<Value>, keep_root: Option<Value>| {
      vec![
        SettingPatch {
          key: o!("reroot"),
          value: reroot.map(JsonValue),
        },
        SettingPatch {
          key: o!("outgroup"),
          value: outgroup.map(JsonValue),
        },
        SettingPatch {
          key: o!("keep_root"),
          value: keep_root.map(JsonValue),
        },
      ]
    };
    let root = command_settings(command)
      .unwrap()
      .choices
      .into_iter()
      .find(|choice| choice.choice == ChoiceName::Root);
    assert_eq!(
      Some(SettingChoice {
        choice: ChoiceName::Root,
        keys: vec_of_owned!["reroot", "outgroup", "keep_root"],
        options: vec![
          ChoiceOption {
            option: ChoiceOptionName::Reroot,
            patch: patch(reroot_method, None, None),
            settings: vec_of_owned!["reroot"],
          },
          ChoiceOption {
            option: ChoiceOptionName::Outgroup,
            patch: patch(None, Some(json!(["outgroup_tip"])), None),
            settings: vec_of_owned!["outgroup"],
          },
          ChoiceOption {
            option: ChoiceOptionName::Keep,
            patch: patch(None, None, Some(json!(true))),
            settings: vec![],
          },
        ],
      }),
      root
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree(       AppCommand::Timetree, json!({}),                     ChoiceOptionName::Reroot)]
  #[case::clock(          AppCommand::Clock,    json!({}),                     ChoiceOptionName::Reroot)]
  #[case::optimize(       AppCommand::Optimize, json!({}),                     ChoiceOptionName::Keep)]
  #[case::optimize_method(AppCommand::Optimize, json!({ "reroot": "min-dev" }), ChoiceOptionName::Reroot)]
  #[trace]
  fn test_choices_root_without_rooting_settings_reports_what_the_command_does(
    #[case] command: AppCommand,
    #[case] rooting: Value,
    #[case] expected: ChoiceOptionName,
  ) {
    let mut draft = fresh_draft(command);
    draft.extend(rooting.as_object().unwrap().clone());

    let root = active_choices_of(command, &Value::Object(draft))
      .into_iter()
      .find(|active| active.choice == ChoiceName::Root);

    assert_eq!(
      Some(ActiveChoice {
        choice: ChoiceName::Root,
        option: expected,
      }),
      root
    );
  }

  #[test]
  fn test_choices_picking_an_option_makes_the_check_report_it() {
    let mismatches = AppCommand::iter()
      .flat_map(|command| {
        command_settings(command)
          .unwrap()
          .choices
          .into_iter()
          .flat_map(move |choice| {
            choice.options.into_iter().filter_map(move |option| {
              let mut draft = fresh_draft(command);
              for patch in &option.patch {
                match &patch.value {
                  Some(value) => {
                    draft.insert(patch.key.clone(), value.0.clone());
                  },
                  None => {
                    draft.remove(&patch.key);
                  },
                }
              }
              let expected = ActiveChoice {
                choice: choice.choice,
                option: option.option,
              };
              let reported = active_choices_of(command, &Value::Object(draft));
              (!reported.contains(&expected)).then(|| format!("{command}: {expected:?} not in {reported:?}"))
            })
          })
      })
      .collect_vec();
    assert_eq!(Vec::<String>::new(), mismatches);
  }

  #[test]
  fn test_choices_invalid_config_that_merges_still_reports_its_options() {
    let response = check(
      AppCommand::Timetree,
      &json!({ "tree": "t.nwk", "metadata": "m.tsv", "keep_root": true, "max_iter": "many" }),
    );
    let CheckConfigResponse::Invalid { choices, .. } = response else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!(
      vec![
        ActiveChoice {
          choice: ChoiceName::ClockRate,
          option: ChoiceOptionName::Estimate
        },
        ActiveChoice {
          choice: ChoiceName::Root,
          option: ChoiceOptionName::Keep
        },
      ],
      choices
    );
  }

  mod helpers {
    use crate::check_config::{CheckConfigRequest, CheckConfigResponse, check_config};
    use crate::config::choices::ActiveChoice;
    use app_commands::command::AppCommand;
    use app_commands::json_value::SparseConfig;
    use serde_json::Value;

    pub(super) fn check(command: AppCommand, config: &Value) -> CheckConfigResponse {
      check_config(&CheckConfigRequest {
        command,
        text: serde_json::to_string(config).unwrap(),
        inputs: SparseConfig::default(),
        input_facts: None,
        folder: None,
        replace_inputs: SparseConfig::default(),
      })
    }

    pub(super) fn active_choices_of(command: AppCommand, config: &Value) -> Vec<ActiveChoice> {
      match check(command, config) {
        CheckConfigResponse::Valid { choices, .. } | CheckConfigResponse::Invalid { choices, .. } => choices,
      }
    }
  }
}
