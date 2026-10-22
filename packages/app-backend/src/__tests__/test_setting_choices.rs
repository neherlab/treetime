#[cfg(test)]
mod tests {
  use crate::__tests__::test_check_config::tests::helpers::fresh_draft;
  use crate::check_config::CheckConfigResponse;
  use crate::config::catalog::command_settings;
  use crate::config::choices::{ActiveChoice, ChoiceName, ChoiceOptionName};
  use app_commands::command::AppCommand;
  use app_commands::config::properties::setting_properties;
  use helpers::{active_choices_of, check};
  use itertools::{Itertools, iproduct};
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use std::collections::BTreeSet;
  use strum::IntoEnumIterator;
  use treetime::timetree::coalescent_timescale::{CoalescentMode, coalescent_mode};

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
        (
          AppCommand::Timetree,
          vec![ChoiceName::ClockRate, ChoiceName::CoalescentPrior, ChoiceName::Root]
        ),
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
  fn test_choices_active_coalescent_prior_follows_the_core_precedence() {
    let reported = iproduct!([None, Some(0.5)], [false, true], [false, true])
      .filter_map(|(coalescent, coalescent_opt, coalescent_skyline)| {
        let mut draft = fresh_draft(AppCommand::Timetree);
        if let Some(tc) = coalescent {
          draft.insert("coalescent".to_owned(), json!(tc));
        }
        draft.insert("coalescent_opt".to_owned(), json!(coalescent_opt));
        draft.insert("coalescent_skyline".to_owned(), json!(coalescent_skyline));
        let response = check(AppCommand::Timetree, &Value::Object(draft));
        let CheckConfigResponse::Valid { choices, .. } = response else {
          return None;
        };
        let option = choices
          .iter()
          .find(|choice| choice.choice == ChoiceName::CoalescentPrior)
          .map(|choice| choice.option);
        let expected = match coalescent_mode(coalescent, coalescent_opt, coalescent_skyline) {
          CoalescentMode::Disabled => ChoiceOptionName::None,
          CoalescentMode::Fixed(_) => ChoiceOptionName::Fixed,
          CoalescentMode::Constant => ChoiceOptionName::Optimized,
          CoalescentMode::Skyline => ChoiceOptionName::Skyline,
        };
        Some((option == Some(expected), expected))
      })
      .collect_vec();
    let accepted_options: BTreeSet<String> = reported.iter().map(|(_, option)| format!("{option:?}")).collect();
    assert_eq!(
      (true, 4),
      (reported.iter().all(|(matches, _)| *matches), accepted_options.len())
    );
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
          choice: ChoiceName::CoalescentPrior,
          option: ChoiceOptionName::None
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
