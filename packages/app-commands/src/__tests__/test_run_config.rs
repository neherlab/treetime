#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use helpers::{ZIKA_METADATA, ZIKA_TREE, resolve};
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::path::Path;

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

  #[test]
  fn test_run_config_hash_depends_on_input_contents_not_paths() {
    let relative = resolve(
      AppCommand::Clock,
      json!({ "tree": ZIKA_TREE, "metadata": ZIKA_METADATA }),
    );
    let absolute = resolve(
      AppCommand::Clock,
      json!({
        "tree": Path::new(ZIKA_TREE).canonicalize().unwrap(),
        "metadata": Path::new(ZIKA_METADATA).canonicalize().unwrap(),
      }),
    );
    let changed = resolve(
      AppCommand::Clock,
      json!({ "tree": ZIKA_TREE, "metadata": ZIKA_METADATA, "keep_root": true }),
    );
    assert_eq!(
      (true, true, false),
      (
        relative["config_hash"].is_string(),
        relative["config_hash"] == absolute["config_hash"],
        relative["config_hash"] == changed["config_hash"]
      )
    );
  }

  #[test]
  fn test_run_config_hash_is_absent_when_an_input_cannot_be_read() {
    let response = resolve(
      AppCommand::Clock,
      json!({ "tree": "missing.nwk", "metadata": ZIKA_METADATA }),
    );
    assert_eq!(
      (
        json!(null),
        json!(
          "When hashing the input of setting `tree`: When opening input 'missing.nwk': No such file or directory (os error 2)"
        )
      ),
      (response["config_hash"].clone(), response["config_hash_error"].clone())
    );
  }

  #[test]
  fn test_run_config_command_line_includes_the_run_outputs() {
    let response = resolve(
      AppCommand::Timetree,
      json!({ "tree": "t.nwk", "metadata": "m.tsv", "output_selection": ["Nwk"] }),
    );
    assert_eq!(
      json!(
        "treetime timetree \\\n  --tree t.nwk \\\n  --metadata m.tsv \\\n  --output-selection 'nwk,auspice,tracelog' \\\n  --output-all out"
      ),
      response["code"]["command_line_text"]
    );
  }

  mod helpers {
    use crate::command::AppCommand;
    use crate::run_config::{RunConfigRequest, run_config};
    use serde_json::Value;

    pub(super) const ZIKA_TREE: &str = "../../data/zika/20/tree.nwk";

    pub(super) const ZIKA_METADATA: &str = "../../data/zika/20/metadata.tsv";

    pub(super) fn resolve(command: AppCommand, config: Value) -> Value {
      let response = run_config(
        &RunConfigRequest { command, config },
        Box::new(|_config: &mut Value| Ok(())),
      );
      serde_json::to_value(response).unwrap()
    }
  }
}
