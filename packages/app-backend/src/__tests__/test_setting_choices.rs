#[cfg(test)]
mod tests {
  use crate::__tests__::test_check_config::tests::helpers::fresh_draft;
  use crate::check_config::CheckConfigResponse;
  use crate::command::offered_commands;
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
  use treetime_utils::{o, vec_of_owned};

  #[test]
  fn test_choices_name_only_settings_of_their_command() {
    let unknown = offered_commands()
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
          .chain(
            settings
              .main_rows
              .iter()
              .flat_map(|row| row.settings.iter().map(|setting| setting.key.clone())),
          )
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
    let choices = offered_commands()
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
            label: o!("Reroot"),
            note: None,
            patch: patch(reroot_method, None, None),
            settings: vec_of_owned!["reroot"],
          },
          ChoiceOption {
            option: ChoiceOptionName::Outgroup,
            label: o!("Outgroup"),
            note: Some(o!("The root is placed on the branch above the most recent common ancestor of the named tips.")),
            patch: patch(None, None, None),
            settings: vec_of_owned!["outgroup"],
          },
          ChoiceOption {
            option: ChoiceOptionName::Keep,
            label: o!("Keep root"),
            note: Some(o!("The input root is kept; use this when the input tree is already rooted.")),
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

    let root = active_choices_of(command, &Value::Object(draft), &[])
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
    let mismatches = offered_commands()
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
              let reported = active_choices_of(command, &Value::Object(draft), &[expected]);
              (!reported.contains(&expected)).then(|| format!("{command}: {expected:?} not in {reported:?}"))
            })
          })
      })
      .collect_vec();
    assert_eq!(Vec::<String>::new(), mismatches);
  }

  #[test]
  fn test_choices_never_fill_a_tip_name() {
    let filled = offered_commands()
      .flat_map(|command| command_settings(command).unwrap().choices)
      .flat_map(|choice| choice.options)
      .flat_map(|option| option.patch)
      .filter(|patch| patch.key == "outgroup" && patch.value.is_some())
      .collect_vec();
    assert_eq!(Vec::<SettingPatch>::new(), filled);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::outgroup_without_tips(ChoiceOptionName::Outgroup, json!({}),                      ChoiceOptionName::Outgroup)]
  #[case::outgroup_with_tips(   ChoiceOptionName::Outgroup, json!({ "outgroup": ["A"] }),   ChoiceOptionName::Outgroup)]
  #[case::reroot_with_tips(     ChoiceOptionName::Reroot,   json!({ "outgroup": ["A"] }),   ChoiceOptionName::Outgroup)]
  #[case::outgroup_kept_root(   ChoiceOptionName::Outgroup, json!({ "keep_root": true }),   ChoiceOptionName::Keep)]
  #[trace]
  fn test_choices_picked_option_counts_while_the_config_agrees(
    #[case] picked: ChoiceOptionName,
    #[case] rooting: Value,
    #[case] expected: ChoiceOptionName,
  ) {
    let mut draft = fresh_draft(AppCommand::Timetree);
    draft.extend(rooting.as_object().unwrap().clone());
    let picked = ActiveChoice { choice: ChoiceName::Root, option: picked };

    let root = active_choices_of(AppCommand::Timetree, &Value::Object(draft), &[picked])
      .into_iter()
      .find(|active| active.choice == ChoiceName::Root);

    assert_eq!(Some(ActiveChoice { choice: ChoiceName::Root, option: expected }), root);
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
      check_picked(command, config, &[])
    }

    pub(super) fn check_picked(command: AppCommand, config: &Value, picked: &[ActiveChoice]) -> CheckConfigResponse {
      check_config(&CheckConfigRequest {
        command,
        text: serde_json::to_string(config).unwrap(),
        inputs: SparseConfig::default(),
        input_facts: None,
        folder: None,
        replace_inputs: SparseConfig::default(),
        choices: picked.to_vec(),
        run_id: None,
      })
    }

    pub(super) fn active_choices_of(command: AppCommand, config: &Value, picked: &[ActiveChoice]) -> Vec<ActiveChoice> {
      match check_picked(command, config, picked) {
        CheckConfigResponse::Valid { choices, .. } | CheckConfigResponse::Invalid { choices, .. } => choices,
      }
    }
  }
}
