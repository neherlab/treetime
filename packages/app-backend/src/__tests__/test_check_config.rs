#[cfg(test)]
pub(crate) mod tests {
  use crate::check_config::{CheckConfigRequest, CheckConfigResponse, check_config};
  use crate::check_inputs::{AlignmentFacts, InputFacts, TreeFacts};
  use crate::run_checks::{CheckLevel, RunCheck};
  use app_commands::command::AppCommand;
  use app_commands::config::source::{ConfigProblem, ConfigSpan};
  use app_commands::json_value::SparseConfig;
  use app_datasets::schema_directive;
  use treetime_utils::io::json::from_json_value;
  use treetime_utils::io::json::to_json_value;

  use helpers::fresh_draft;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Map, json};
  use std::path::PathBuf;
  use treetime_utils::assert_error;
  use treetime_utils::o;

  const SEQUENCE_LENGTH_REQUIRED: &str =
    "--sequence-length is required without an alignment: the time inference needs the sequence length";

  const SEQUENCE_LENGTH_WITH_ALIGNMENT: &str = "--sequence-length applies only without --alignment: the run takes \
                                                the sequence length from the alignment, and an alignment of \
                                                variable sites with a separate genome length is not supported";

  const MARGINAL_NEEDS_ALIGNMENT: &str =
    "--branch-length-mode=marginal requires an alignment: add --alignment or use --branch-length-mode=input";

  #[test]
  fn test_check_config_valid_yaml_returns_config_with_defaults() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Timetree,
      text: indoc! {r#"
        "$schema": "https://example.org/input-config-timetree.schema.json"
        tree: data/zika/20/tree.nwk
        metadata: data/zika/20/metadata.tsv
        alignment: [data/zika/20/aln.fasta.xz]
        max_iter: 4
      "#}
      .to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig::default(),
    });
    let CheckConfigResponse::Valid { config, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(
      (json!("data/zika/20/tree.nwk"), json!(4), json!(4.0), None),
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
        metadata: data/zika/20/metadata.tsv
        alignment: [data/zika/20/aln.fasta.xz]
        output_all: /home/user/results
        output_tree_nwk: /home/user/tree.nwk
      "#}
      .to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig::default(),
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
      text: r#"{ "tree": "t.nwk", "alignment": ["a.fasta"], "method_anc": "parsimony" }"#.to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig::default(),
    });
    let CheckConfigResponse::Valid { config, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(json!("parsimony"), config["method_anc"]);
  }

  #[test]
  fn test_check_config_a_setting_without_effect_blocks_at_both_settings() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Ancestral,
      text: "tree: t.nwk\nalignment: [a.fasta]\nmethod_anc: parsimony\nmodel: infer\n".to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig::default(),
    });
    let CheckConfigResponse::Invalid { checks, .. } = response else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!(
      vec![(
        o!("config-0"),
        o!("--model applies only with --method-anc marginal: parsimony reconstruction uses no substitution model"),
        vec![o!("model"), o!("method_anc")]
      )],
      checks
        .into_iter()
        .map(|check| (check.id, check.text, check.settings))
        .collect::<Vec<_>>()
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::ancestral_parsimony( AppCommand::Ancestral, json!({ "method_anc": "parsimony" }))]
  #[case::homoplasy_parsimony( AppCommand::Homoplasy, json!({ "method_anc": "parsimony" }))]
  #[case::optimize_dense(      AppCommand::Optimize,  json!({ "representation": "dense" }))]
  #[case::timetree_opt(        AppCommand::Timetree,  json!({ "coalescent": "opt" }))]
  #[case::timetree_input_mode( AppCommand::Timetree,  json!({ "branch_length_mode": "input" }))]
  #[trace]
  fn test_check_config_the_defaults_of_a_fresh_draft_give_no_setting(
    #[case] command: AppCommand,
    #[case] settings: serde_json::Value,
  ) {
    let mut draft = fresh_draft(command);
    draft.extend(settings.as_object().unwrap().clone());
    let response = check_config(&CheckConfigRequest {
      command,
      text: serde_json::Value::Object(draft).to_string(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig::default(),
    });
    assert!(
      matches!(response, CheckConfigResponse::Valid { .. }),
      "expected a valid config, got {response:?}"
    );
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
      replace_inputs: SparseConfig::default(),
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
      replace_inputs: SparseConfig::default(),
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
      replace_inputs: SparseConfig::default(),
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
      text: "tree: tree.nwk\nalignment: [a.fasta]\n".to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig::default(),
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
      AppCommand::Homoplasy,
      AppCommand::Clock,
      AppCommand::Mugration
    )]
    command: AppCommand,
  ) {
    let draft = fresh_draft(command);
    let response = check_config(&CheckConfigRequest {
      command,
      text: serde_json::to_string(&draft).unwrap(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig::default(),
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
      replace_inputs: SparseConfig::default(),
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
      replace_inputs: SparseConfig::default(),
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
        "the following required arguments were not provided:\n  --metadata <FILE>\n  --attribute <COLUMN>",
        0,
        None
      ),
      (message.as_str(), problems.len(), rendered)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::alignment(                   "alignment: [a.fasta]",                                None)]
  #[case::auto_sequence_length(        "sequence_length: 1000",                               None)]
  #[case::input_sequence_length(       "branch_length_mode: input\nsequence_length: 1000",     None)]
  #[case::auto_neither(                "max_iter: 2",                                         Some((["alignment", "sequence_length"],    SEQUENCE_LENGTH_REQUIRED)))]
  #[case::input_neither(               "branch_length_mode: input",                           Some((["alignment", "sequence_length"],    SEQUENCE_LENGTH_REQUIRED)))]
  #[case::zero_sequence_length(        "branch_length_mode: input\nsequence_length: 0",        Some((["alignment", "sequence_length"],    SEQUENCE_LENGTH_REQUIRED)))]
  #[case::marginal_without_alignment(  "branch_length_mode: marginal\nsequence_length: 1000",  Some((["alignment", "branch_length_mode"], MARGINAL_NEEDS_ALIGNMENT)))]
  #[case::alignment_and_length(        "alignment: [a.fasta]\nsequence_length: 1000",          Some((["sequence_length", "alignment"],    SEQUENCE_LENGTH_WITH_ALIGNMENT)))]
  #[trace]
  fn test_check_config_timetree_blocks_a_run_without_sequence_length(
    #[case] settings: &str,
    #[case] expected: Option<([&str; 2], &str)>,
  ) {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Timetree,
      text: format!("tree: t.nwk\nmetadata: m.tsv\n{settings}\n"),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig::default(),
    });
    let expected = expected
      .map(|(settings, text)| RunCheck {
        id: o!("config-0"),
        level: CheckLevel::Block,
        text: o!(text),
        settings: settings.iter().map(|setting| (*setting).to_owned()).collect(),
        fix: None,
      })
      .into_iter()
      .collect::<Vec<_>>();
    let checks = match response {
      CheckConfigResponse::Valid { checks, .. } | CheckConfigResponse::Invalid { checks, .. } => checks,
    };
    assert_eq!(expected, checks);
  }

  #[test]
  fn test_check_config_yaml_syntax_error_is_invalid() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Clock,
      text: "tree: t.nwk\ntree: other.nwk\n".to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig::default(),
    });
    let CheckConfigResponse::Invalid { message, problems, .. } = response else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!(
      (
        "invalid configuration: could not parse config: Unexpected: duplicate map key \"tree\" at line 2 column 1 (path: tree)",
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
      replace_inputs: SparseConfig::default(),
    });
    let value = to_json_value(&response).unwrap();
    assert_eq!(
      (json!("valid"), json!("prune")),
      (value["status"].clone(), value["command"].clone())
    );
  }

  #[test]
  fn test_check_config_request_rejects_unknown_command() {
    let result = from_json_value::<CheckConfigRequest>(&json!({ "command": "transmission", "text": "" }));
    assert_error!(
      result,
      "When converting a JSON value: Unexpected: unknown variant `transmission` of AppCommand, expected one of `timetree`, `clock`, `ancestral`, `homoplasy`, `mugration`, `optimize`, `prune`"
    );
  }

  #[test]
  fn test_check_config_schema_directive_selects_the_command() {
    let text = format!(
      "{}\ntree: t.nwk\nalignment: [a.fasta]\ncollapse_empty: true\n",
      schema_directive("prune")
    );
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Timetree,
      text,
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig::default(),
    });
    let CheckConfigResponse::Valid { command, config, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(
      (AppCommand::Prune, json!(true)),
      (command, config["collapse_empty"].clone())
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
      replace_inputs: SparseConfig::default(),
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
      text: "alignment: [\"text.fasta\"]\nrepresentation: dense".to_owned(),
      inputs: SparseConfig(inputs.as_object().unwrap().clone()),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig::default(),
    });
    let CheckConfigResponse::Valid { config, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(
      (json!("draft.nwk"), json!(["text.fasta"]), json!("dense")),
      (
        config["tree"].clone(),
        config["alignment"].clone(),
        config["representation"].clone()
      )
    );
  }

  #[test]
  fn test_check_config_replace_inputs_swaps_the_paths_of_an_example_config() {
    let replacements = json!({
      "tree": "inputs/1/tree.nwk",
      "metadata": "inputs/2/metadata.tsv",
      "weights": "inputs/3/country_weights.csv",
    });
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Mugration,
      text: include_str!("../../../../data/zika/20/mugration-country-weights.yaml").to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig(replacements.as_object().unwrap().clone()),
    });
    let CheckConfigResponse::Valid { config, code, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(
      (
        json!("inputs/1/tree.nwk"),
        json!("inputs/2/metadata.tsv"),
        json!("inputs/3/country_weights.csv"),
        json!("country"),
        true
      ),
      (
        config["tree"].clone(),
        config["metadata"].clone(),
        config["weights"].clone(),
        config["attribute"].clone(),
        code
          .command_line
          .iter()
          .any(|line| line.text == "--weights inputs/3/country_weights.csv")
      )
    );
  }

  #[test]
  fn test_check_config_replace_inputs_rejects_a_setting_that_is_not_an_input() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Mugration,
      text: "tree: tree.nwk\nmetadata: metadata.tsv\nattribute: country\n".to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig(json!({ "attribute": "region" }).as_object().unwrap().clone()),
    });
    let CheckConfigResponse::Invalid { message, .. } = response else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!(
      "`replace_inputs` names `attribute`, which is not an input setting of the command",
      message
    );
  }

  #[test]
  fn test_check_config_replace_inputs_checks_the_type_of_the_values() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Mugration,
      text: "tree: tree.nwk\nmetadata: metadata.tsv\nattribute: country\n".to_owned(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig(json!({ "tree": 5 }).as_object().unwrap().clone()),
    });
    let CheckConfigResponse::Invalid { problems, .. } = response else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!(
      vec![o!("config::type")],
      problems.into_iter().map(|problem| problem.code).collect::<Vec<_>>()
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
      replace_inputs: SparseConfig::default(),
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
      replace_inputs: SparseConfig::default(),
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
      alignments: vec![AlignmentFacts {
        number: 1,
        path: PathBuf::from("a.fasta"),
        sequences: 2,
        min_length: 4,
        max_length: 4,
        duplicate_names: vec![],
        tips_without_sequence: Some(vec![o!("A")]),
      }],
      tree: Some(TreeFacts {
        tips: 3,
        internal_nodes: 2,
        polytomies: 0,
        unnamed_tips: 0,
        duplicate_node_names: vec![],
      }),
      ..InputFacts::default()
    };
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Ancestral,
      text: "tree: t.nwk\nalignment: [a.fasta]\n".to_owned(),
      inputs: SparseConfig::default(),
      input_facts: Some(facts),
      folder: None,
      replace_inputs: SparseConfig::default(),
    });
    let CheckConfigResponse::Valid { checks, code, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(
      (
        vec![o!(
          "1 of 3 tree tips have no sequence in the alignment: A. TreeTime treats their sequences as unknown."
        )],
        o!("treetime ancestral \\\n  --alignment a.fasta \\\n  --tree t.nwk \\\n  --output-all out")
      ),
      (
        checks.into_iter().map(|check| check.text).collect::<Vec<_>>(),
        code.command_line_text
      )
    );
  }

  pub(crate) mod helpers {
    use crate::check_inputs::{InputKind, InputNeed};
    use crate::command::AppCommandExt;
    use crate::config::catalog::command_settings;
    use app_commands::command::AppCommand;
    use serde_json::{Map, Value, json};
    use treetime_utils::o;

    pub(crate) fn fresh_draft(command: AppCommand) -> Map<String, Value> {
      let mut draft = Map::new();
      for spec in command_settings(command).unwrap().settings {
        if let Some(default) = spec.default_value {
          draft.insert(spec.key, default.0);
        }
      }
      for input in command
        .inputs()
        .iter()
        .filter(|input| input.need != InputNeed::Optional)
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
      draft
    }
  }
}
