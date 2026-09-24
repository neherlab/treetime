#[cfg(test)]
mod tests {
  use crate::command::{AppCommand, CheckConfigRequest, CheckConfigResponse, check_config};
  use crate::config::source::{ConfigProblem, ConfigSpan};
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use treetime_utils::assert_error;

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
    });
    let CheckConfigResponse::Valid { config } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    assert_eq!(
      (json!("data/zika/20/tree.nwk"), json!(4), json!(3.0), None),
      (
        config["tree"].clone(),
        config["max_iter"].clone(),
        config["clock_filter"].clone(),
        config.get("$schema")
      )
    );
  }

  #[test]
  fn test_check_config_accepts_json_text() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Ancestral,
      text: r#"{ "tree": "t.nwk", "method_anc": "parsimony" }"#.to_owned(),
    });
    let CheckConfigResponse::Valid { config } = response else {
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
    });
    let CheckConfigResponse::Invalid {
      message,
      causes,
      problems,
      rendered,
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

  #[test]
  fn test_check_config_wrong_enum_value_is_invalid() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Ancestral,
      text: "tree: t.nwk\nmethod_anc: margnal\n".to_owned(),
    });
    let CheckConfigResponse::Invalid { message, problems, .. } = response else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!(
      (
        "invalid configuration: `margnal` is not a valid value",
        vec!["config::enum".to_owned()],
        vec![Some(
          "did you mean `marginal`? Valid values: `joint`, `marginal`, `parsimony`".to_owned()
        )],
      ),
      (
        message.as_str(),
        problems.iter().map(|problem| problem.code.clone()).collect::<Vec<_>>(),
        problems.iter().map(|problem| problem.help.clone()).collect::<Vec<_>>(),
      )
    );
  }

  #[test]
  fn test_check_config_missing_required_input_uses_cli_wording() {
    let response = check_config(&CheckConfigRequest {
      command: AppCommand::Mugration,
      text: "tree: t.nwk\n".to_owned(),
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
    let response = CheckConfigResponse::Valid { config: json!({}) };
    assert_eq!(
      json!({ "status": "valid", "config": {} }),
      serde_json::to_value(response).unwrap()
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
}
