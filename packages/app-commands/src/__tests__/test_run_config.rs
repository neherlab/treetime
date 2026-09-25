#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use helpers::resolve;
  use pretty_assertions::assert_eq;
  use serde_json::json;

  #[test]
  fn test_run_config_adds_the_run_outputs_to_the_chosen_selection() {
    let response = resolve(
      AppCommand::Timetree,
      json!({ "tree": "t.nwk", "metadata": "m.tsv", "output_selection": ["Nwk"], "output_all": "/elsewhere" }),
    );
    assert_eq!(
      (json!("valid"), json!(["Nwk", "Auspice", "Tracelog"]), json!("out")),
      (
        response["status"].clone(),
        response["config"]["output_selection"].clone(),
        response["config"]["output_all"].clone()
      )
    );
  }

  #[test]
  fn test_run_config_reports_an_unknown_setting_as_the_cli_does() {
    let response = resolve(
      AppCommand::Timetree,
      json!({ "tree": "t.nwk", "metadata": "m.tsv", "prune_short": 1 }),
    );
    assert_eq!(
      (
        json!("invalid"),
        json!("invalid configuration: unknown field `prune_short`")
      ),
      (response["status"].clone(), response["message"].clone())
    );
  }

  mod helpers {
    use crate::command::AppCommand;
    use crate::run_config::{RunConfigRequest, run_config};
    use serde_json::Value;

    pub(super) fn resolve(command: AppCommand, config: Value) -> Value {
      serde_json::to_value(run_config(&RunConfigRequest { command, config })).unwrap()
    }
  }
}
