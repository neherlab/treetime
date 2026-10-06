#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use helpers::{Source, apply};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};

  const SETTING_5000: &str = "analysis:\n  max_grid_points: 5000\n";

  const NO_SETTING: &str = "ui:\n  theme: dark\n";

  const BROKEN: &str = "analysis: [\n";

  #[rustfmt::skip]
  #[rstest]
  #[case::unset_leaves_config(          Source::Unset,                  AppCommand::Timetree, json!({}),                          json!({}))]
  #[case::settings_fill_unset_value(    Source::Settings(SETTING_5000), AppCommand::Timetree, json!({}),                          json!({ "max_grid_points": 5000 }))]
  #[case::settings_keep_lower_value(    Source::Settings(SETTING_5000), AppCommand::Timetree, json!({ "max_grid_points": 3000 }), json!({ "max_grid_points": 3000 }))]
  #[case::settings_keep_higher_value(   Source::Settings(SETTING_5000), AppCommand::Timetree, json!({ "max_grid_points": 9000 }), json!({ "max_grid_points": 9000 }))]
  #[case::settings_without_the_setting( Source::Settings(NO_SETTING),   AppCommand::Timetree, json!({}),                          json!({}))]
  #[case::settings_not_read_for_clock(  Source::Settings(BROKEN),       AppCommand::Clock,    json!({}),                          json!({}))]
  #[case::settings_not_read_when_set(   Source::Settings(BROKEN),       AppCommand::Timetree, json!({ "max_grid_points": 3000 }), json!({ "max_grid_points": 3000 }))]
  #[case::server_fills_unset_value(     Source::Server(4000),           AppCommand::Timetree, json!({}),                          json!({ "max_grid_points": 4000 }))]
  #[case::server_keeps_lower_value(     Source::Server(4000),           AppCommand::Timetree, json!({ "max_grid_points": 2000 }), json!({ "max_grid_points": 2000 }))]
  #[case::server_keeps_equal_value(     Source::Server(4000),           AppCommand::Timetree, json!({ "max_grid_points": 4000 }), json!({ "max_grid_points": 4000 }))]
  #[case::server_skips_clock(           Source::Server(4000),           AppCommand::Clock,    json!({}),                          json!({}))]
  #[trace]
  fn test_run_limits_apply(
    #[case] source: Source,
    #[case] command: AppCommand,
    #[case] config: Value,
    #[case] expected: Value,
  ) {
    assert_eq!(Ok(expected), apply(source, command, config));
  }

  #[test]
  fn test_run_limits_server_rejects_a_higher_value() {
    let expected = "`max_grid_points` is 4001, more than the limit of this server, 4000 (treetime-server \
      --max-grid-points); set 4000 or less, or leave the setting unset";
    assert_eq!(
      Err(expected.to_owned()),
      apply(
        Source::Server(4000),
        AppCommand::Timetree,
        json!({ "max_grid_points": 4001 })
      )
    );
  }

  #[test]
  fn test_run_limits_unreadable_settings_file_fails() {
    let actual = apply(Source::Settings(BROKEN), AppCommand::Timetree, json!({}));
    assert!(
      actual
        .as_ref()
        .is_err_and(|message| message.starts_with("When reading the settings file '")),
      "{actual:?}"
    );
  }

  mod helpers {
    use crate::app_settings::store::{AppSettingsStore, SETTINGS_YAML};
    use crate::command::AppCommand;
    use crate::run_limits::RunLimits;
    use serde_json::Value;
    use std::fs;
    use std::sync::Arc;
    use tempfile::tempdir;
    use treetime_grid::MaxGridPoints;
    use treetime_utils::error::report_to_string;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Source {
      Unset,
      Settings(&'static str),
      Server(usize),
    }

    pub(super) fn apply(source: Source, command: AppCommand, config: Value) -> Result<Value, String> {
      let dir = tempdir().unwrap();
      let limits = match source {
        Source::Unset => RunLimits::Unset,
        Source::Settings(text) => {
          fs::write(dir.path().join(SETTINGS_YAML), text).unwrap();
          RunLimits::Settings(Arc::new(AppSettingsStore::open(dir.path()).unwrap()))
        },
        Source::Server(limit) => RunLimits::Server(MaxGridPoints::new(limit).unwrap()),
      };
      let Value::Object(mut settings) = config else {
        panic!("a test config must be a mapping");
      };
      limits
        .apply(command, &mut settings)
        .map_err(|report| report_to_string(&report))?;
      Ok(Value::Object(settings))
    }
  }
}
