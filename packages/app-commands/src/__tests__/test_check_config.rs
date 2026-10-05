#[cfg(test)]
mod tests {
  use crate::check_config::{CheckConfigRequest, CheckConfigResponse, check_config};
  use crate::check_inputs::{InputFacts, InputKind, InputNeed, TreeFacts};
  use crate::command::AppCommand;
  use crate::config::catalog::command_settings;
  use crate::config::source::{ConfigProblem, ConfigSpan};
  use crate::json_value::SparseConfig;
  use crate::run_checks::CheckLevel;
  use app_datasets::schema_directive;
  use eyre::Report;
  use helpers::set_at;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Map, json};
  use std::path::PathBuf;
  use treetime_utils::assert_error;
  use treetime_utils::o;

  #[test]
  fn test_check_config_valid_yaml_returns_config_with_defaults() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Timetree,
      text: indoc! {r#"
        "$schema": "https://example.org/input-config-timetree.schema.json"
        tree: data/zika/20/tree.nwk
        max_iter: 4
      "#}
      .to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Valid { config, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(
      (json!("data/zika/20/tree.nwk"), json!(4), json!(3.0), None),
      (
        config["tree"].clone(),
        config["max_iter"].clone(),
        config["clock_filter"].clone(),
        config.get("$schema").cloned()
      )
    );
  }

  #[test]
  fn test_check_config_drops_the_output_paths_the_app_sets_for_each_run() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Timetree,
      text: indoc! {r#"
        tree: data/zika/20/tree.nwk
        output_all: /home/user/results
        output_tree_nwk: /home/user/tree.nwk
      "#}
      .to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Valid { config, code, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(
      (None, None, Some("--output-all out")),
      (
        config.get("output_all").cloned(),
        config.get("output_tree_nwk").cloned(),
        code.command_line.last().map(|line| line.text.as_str())
      )
    );
  }

  #[test]
  fn test_check_config_accepts_json_text() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Ancestral,
      text: r#"{ "tree": "t.nwk", "method_anc": "parsimony" }"#.to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Valid { config, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(json!("parsimony"), config["method_anc"]);
  }

  #[test]
  fn test_check_config_unknown_field_is_invalid_with_location() {
    let text = "tree: t.nwk\ndefinitely_not_a_real_field: 1\n";
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Ancestral,
      text: text.to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Invalid {
      message,
      causes,
      problems,
      rendered,
      ..
    } = response
    else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!(
      (
        "invalid configuration: unknown field `definitely_not_a_real_field`",
        vec![],
        vec![ConfigProblem {
          code: "config::unknown-field".to_owned(),
          message: "unknown field `definitely_not_a_real_field`".to_owned(),
          span: Some(ConfigSpan {
            offset: text.find("definitely").unwrap(),
            length: "definitely_not_a_real_field".len(),
          }),
          help: None,
        }],
      ),
      (message.as_str(), causes, problems)
    );
    assert!(rendered.is_some_and(|rendered| rendered.contains("definitely_not_a_real_field: 1")));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::misspelled(   "margnal", "did you mean `marginal`? Valid values: `marginal`, `parsimony`")]
  #[case::removed_joint("joint",   "valid values: `marginal`, `parsimony`")]
  #[trace]
  fn test_check_config_wrong_enum_value_is_invalid(#[case] value: &str, #[case] expected_help: &str) {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Ancestral,
      text: format!("tree: t.nwk\nmethod_anc: {value}\n"),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Invalid { message, problems, .. } = response else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!(
      (
        format!("invalid configuration: `{value}` is not a valid value"),
        vec!["config::enum".to_owned()],
        vec![Some(expected_help.to_owned())],
      ),
      (
        message,
        problems.iter().map(|problem| problem.code.clone()).collect::<Vec<_>>(),
        problems.iter().map(|problem| problem.help.clone()).collect::<Vec<_>>(),
      )
    );
  }

  #[test]
  fn test_check_config_folder_resolves_the_paths_of_the_text_and_keeps_the_draft_inputs() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Ancestral,
      text: "tree: tree.nwk\noutput_all: out\n".to_owned(),
      inputs: SparseConfig(Map::from_iter([(o!("alignment"), json!(["inputs/aln.fasta"]))])),
      input_facts: None,
      folder: Some(PathBuf::from("/data/zika/20")),
    });
    let CheckConfigResponse::Valid { config, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(
      (json!("/data/zika/20/tree.nwk"), json!(["inputs/aln.fasta"]), None),
      (
        config["tree"].clone(),
        config["alignment"].clone(),
        config.get("output_all")
      )
    );
  }

  #[test]
  fn test_check_config_without_folder_keeps_relative_paths() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Ancestral,
      text: "tree: tree.nwk\n".to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Valid { config, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(json!("tree.nwk"), config["tree"]);
  }

  #[rstest]
  #[trace]
  fn test_check_config_accepts_the_fresh_draft_of_each_command(
    #[values(
      AppCommand::Timetree,
      AppCommand::Optimize,
      AppCommand::Prune,
      AppCommand::Ancestral,
      AppCommand::Clock,
      AppCommand::Mugration
    )]
    command: AppCommand,
  ) {
    let mut draft = Map::new();
    for spec in command_settings(command).unwrap().settings {
      if let Some(default) = spec.default_value {
        set_at(&mut draft, &spec.path, default.0);
      }
    }
    for input in command
      .inputs()
      .iter()
      .filter(|input| input.need == InputNeed::Required)
    {
      let (key, value) = match input.kind {
        InputKind::Tree => ("tree", json!("t.nwk")),
        InputKind::Metadata => ("metadata", json!("m.tsv")),
        InputKind::Alignment => ("alignment", json!(["a.fasta"])),
      };
      draft.insert(o!(key), value);
    }
    if command == AppCommand::Mugration {
      draft.insert(o!("attribute"), json!("country"));
    }
    let response = check_config(&CheckConfigRequest {
      command,
      text: serde_json::to_string(&draft).unwrap(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    assert!(
      matches!(response, CheckConfigResponse::Valid { .. }),
      "{command}: {response:?}"
    );
  }

  #[test]
  fn test_check_config_rejects_a_null_value() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Timetree,
      text: "tree: t.nwk\nclock_rate: null\n".to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Invalid { message, .. } = response else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!("invalid configuration: null is not of type \"number\"", message);
  }

  #[test]
  fn test_check_config_missing_required_input_uses_cli_wording() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Mugration,
      text: "tree: t.nwk\n".to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Invalid {
      message,
      problems,
      rendered,
      ..
    } = response
    else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!(
      (
        "the following required arguments were not provided:\n  --metadata <METADATA>\n  --attribute <ATTRIBUTE>",
        0,
        None
      ),
      (message.as_str(), problems.len(), rendered)
    );
  }

  #[test]
  fn test_check_config_yaml_syntax_error_is_invalid() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Clock,
      text: "tree: t.nwk\ntree: other.nwk\n".to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Invalid { message, problems, .. } = response else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!(
      (
        "invalid configuration: could not parse config: error: line 2 column 1: duplicate mapping key: tree, set DuplicateKeyPolicy in Options if acceptable",
        vec!["config::syntax".to_owned()]
      ),
      (
        message.lines().next().unwrap_or_default(),
        problems.iter().map(|problem| problem.code.clone()).collect::<Vec<_>>()
      )
    );
  }

  #[test]
  fn test_check_config_response_serializes_with_status_tag() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Prune,
      text: "tree: t.nwk\n".to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    let value = serde_json::to_value(response).unwrap();
    assert_eq!(
      (json!("valid"), json!("prune")),
      (value["status"].clone(), value["command"].clone())
    );
  }

  #[test]
  fn test_check_config_request_rejects_unknown_command() {
    let result =
      serde_json::from_value::<CheckConfigRequest>(json!({ "command": "homoplasy", "text": "" })).map_err(Report::from);
    assert_error!(
      result,
      "unknown variant `homoplasy`, expected one of `timetree`, `optimize`, `prune`, `ancestral`, `clock`, `mugration`"
    );
  }

  #[test]
  fn test_check_config_schema_directive_selects_the_command() {
    let text = format!("{}\ntree: t.nwk\nprune_empty: true\n", schema_directive("prune"));
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Timetree,
      text,
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Valid { command, config, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(
      (AppCommand::Prune, json!(true)),
      (command, config["prune_empty"].clone())
    );
  }

  #[test]
  fn test_check_config_directive_of_a_command_the_app_does_not_run_keeps_the_requested_command() {
    let text = format!("{}\ntree: t.nwk\n", schema_directive("pipeline"));
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Prune,
      text,
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Valid { command, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(AppCommand::Prune, command);
  }

  #[test]
  fn test_check_config_adds_the_inputs_the_text_does_not_set() {
    let inputs = json!({ "tree": "draft.nwk", "alignment": ["draft.fasta"], "metadata": "draft.tsv" });
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Ancestral,
      text: "alignment: [\"text.fasta\"]\ndense: true".to_owned(),
      inputs: SparseConfig(inputs.as_object().unwrap().clone()),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Valid { config, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(
      (json!("draft.nwk"), json!(["text.fasta"]), json!(true)),
      (
        config["tree"].clone(),
        config["alignment"].clone(),
        config["dense"].clone()
      )
    );
  }

  #[test]
  fn test_check_config_problems_keep_their_places_in_the_text_when_inputs_are_added() {
    let text = "definitely_not_a_real_field: 1";
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Prune,
      text: text.to_owned(),
      inputs: SparseConfig(json!({ "tree": "draft.nwk" }).as_object().unwrap().clone()),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Invalid { problems, .. } = response else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!(
      vec![Some(ConfigSpan {
        offset: 0,
        length: "definitely_not_a_real_field".len(),
      })],
      problems.into_iter().map(|problem| problem.span).collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_check_config_invalid_lists_messages_and_blocking_checks() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Ancestral,
      text: "tree: t.nwk\nmethod_anc: margnal\n".to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
    });
    let CheckConfigResponse::Invalid { messages, checks, .. } = response else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!(
      (
        messages.clone(),
        vec![
          (o!("missing-alignment"), CheckLevel::Block, o!("Add an alignment file.")),
          (o!("config-0"), CheckLevel::Block, messages[0].clone()),
        ]
      ),
      (
        messages,
        checks
          .into_iter()
          .map(|check| (check.id, check.level, check.text))
          .collect::<Vec<_>>()
      )
    );
  }

  #[test]
  fn test_check_config_valid_classifies_the_input_facts() {
    let facts = InputFacts {
      tips_without_sequence: Some(vec![o!("A")]),
      tree: Some(TreeFacts {
        tips: 3,
        internal_nodes: 2,
        polytomies: 0,
        unnamed_tips: 0,
        duplicate_tip_names: vec![],
      }),
      ..InputFacts::default()
    };
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Ancestral,
      text: "tree: t.nwk\nalignment: [a.fasta]\n".to_owned(),
      inputs: SparseConfig::default(),
      input_facts: Some(facts),
      folder: None,
    });
    let CheckConfigResponse::Valid { checks, code, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(
      (
        vec![o!("1 of 3 tree tips have no sequence in the alignment: A")],
        o!("treetime ancestral \\\n  --alignment a.fasta \\\n  --tree t.nwk \\\n  --output-all out")
      ),
      (
        checks.into_iter().map(|check| check.text).collect::<Vec<_>>(),
        code.command_line_text
      )
    );
  }

  mod helpers {
    use serde_json::{Map, Value, json};

    pub(super) fn set_at(config: &mut Map<String, Value>, path: &[String], value: Value) {
      let Some((last, parents)) = path.split_last() else {
        return;
      };
      let parent = parents.iter().fold(config, |map, key| {
        map
          .entry(key.clone())
          .or_insert_with(|| json!({}))
          .as_object_mut()
          .expect("a setting group is a mapping")
      });
      parent.insert(last.clone(), value);
    }
  }
}
