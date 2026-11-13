#[cfg(test)]
mod tests {
  use crate::__tests__::test_check_config::tests::helpers::fresh_draft;
  use crate::check_config::{CheckConfigRequest, CheckConfigResponse, check_config};
  use crate::config::effective::EffectiveValue;
  use app_commands::command::AppCommand;
  use app_commands::json_value::SparseConfig;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree_default(  AppCommand::Timetree, json!({}),                       "Automatic: least-squares rerooting")]
  #[case::clock_default(     AppCommand::Clock,    json!({}),                       "Automatic: least-squares rerooting")]
  #[case::optimize_default(  AppCommand::Optimize, json!({}),                       "Not used: the input root is kept")]
  #[case::outgroup(          AppCommand::Timetree, json!({ "outgroup": ["A"] }),    "Not used: the tree is rooted on the outgroup")]
  #[case::keep_root(         AppCommand::Timetree, json!({ "keep_root": true }),    "Not used: the input root is kept")]
  #[trace]
  fn test_effective_values_rooting(#[case] command: AppCommand, #[case] rooting: Value, #[case] expected: &str) {
    assert_eq!(Some(expected.to_owned()), effective_text(command, &rooting, "reroot"));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::default(   json!({}),                      Some("4 interquartile ranges"))]
  #[case::changed(   json!({ "clock_filter": 3.0 }), None)]
  #[trace]
  fn test_effective_values_clock_filter(#[case] settings: Value, #[case] expected: Option<&str>) {
    assert_eq!(expected.map(str::to_owned), effective_text(AppCommand::Timetree, &settings, "clock_filter"));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::clock_rate(  "clock_rate",   "Estimated from the data")]
  #[case::coalescent(  "coalescent",   "None: the node times have no coalescent prior")]
  #[case::switch(      "confidence",   "Off")]
  #[case::documented(  "n_skyline",    "20")]
  #[case::enum_value(  "model",        "infer")]
  #[case::no_default(  "date_column",  "Not set")]
  #[trace]
  fn test_effective_values_timetree(#[case] key: &str, #[case] expected: &str) {
    assert_eq!(Some(expected.to_owned()), effective_text(AppCommand::Timetree, &json!({}), key));
  }

  #[test]
  fn test_effective_values_leave_out_paths() {
    let keys = effective(AppCommand::Timetree, &json!({}))
      .into_iter()
      .filter(|(key, _)| ["tree", "metadata", "alignment", "output_tree"].contains(&key.as_str()))
      .collect::<Vec<_>>();
    assert_eq!(Vec::<(String, String)>::new(), keys);
  }

  #[test]
  fn test_effective_values_without_inputs() {
    let response = check_config(&request(AppCommand::Timetree, &json!({})));
    let CheckConfigResponse::Invalid { effective, .. } = response else {
      panic!("expected a config without inputs to be invalid, got {response:?}");
    };
    let reroot = effective
      .into_iter()
      .find(|value| value.key == "reroot")
      .map(|value| value.text);
    assert_eq!(Some(o!("Automatic: least-squares rerooting")), reroot);
  }

  #[test]
  fn test_effective_values_of_unparsed_text_are_empty() {
    let response = check_config(&CheckConfigRequest {
      text: o!("tree: [unclosed"),
      ..request(AppCommand::Timetree, &json!({}))
    });
    let CheckConfigResponse::Invalid { effective, .. } = response else {
      panic!("expected an invalid config, got {response:?}");
    };
    assert_eq!(Vec::<EffectiveValue>::new(), effective);
  }

  fn effective_text(command: AppCommand, settings: &Value, key: &str) -> Option<String> {
    effective(command, settings)
      .into_iter()
      .find(|(found, _)| found == key)
      .map(|(_, text)| text)
  }

  fn effective(command: AppCommand, settings: &Value) -> Vec<(String, String)> {
    let mut draft = fresh_draft(command);
    draft.extend(settings.as_object().unwrap().clone());
    let response = check_config(&request(command, &Value::Object(draft)));
    let CheckConfigResponse::Valid { effective, .. } = response else {
      panic!("expected a valid config, got {response:?}");
    };
    effective.into_iter().map(|value| (value.key, value.text)).collect()
  }

  fn request(command: AppCommand, config: &Value) -> CheckConfigRequest {
    CheckConfigRequest {
      command,
      text: serde_json::to_string(config).unwrap(),
      inputs: SparseConfig::default(),
      input_facts: None,
      folder: None,
      replace_inputs: SparseConfig::default(),
      choices: vec![],
      run_id: None,
      example: None,
    }
  }
}
