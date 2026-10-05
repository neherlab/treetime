#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use crate::config::catalog::{SettingRole, command_settings};
  use crate::config::code::{CodeLineKind, config_code};
  use crate::config::settings::setting_mut;
  use helpers::{changed_value, load_yaml, parse_command_line, set_setting, with_output_dir};
  use indoc::indoc;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use serde_json::{Map, Value, json};
  use strum::IntoEnumIterator;

  #[test]
  fn test_config_code_command_line_of_each_setting_parses_back_to_the_same_config() {
    for command in AppCommand::iter() {
      let defaults = command.default_config().unwrap();
      for spec in command_settings(command).unwrap().settings {
        if spec.role == SettingRole::Output {
          continue;
        }
        let mut config = defaults.clone();
        let value = changed_value(&spec);
        set_setting(&mut config, &spec.path, value.clone());
        let code = config_code(command, &config).unwrap();
        assert!(
          code.command_line.iter().all(|line| line.kind != CodeLineKind::Comment),
          "{command}: `{}` = {value} has no command-line form",
          spec.key
        );
        assert_eq!(
          Value::Object(with_output_dir(config)),
          Value::Object(parse_command_line(command, &code.command_line_text)),
          "{command}: `{}` = {value}\n{}",
          spec.key,
          code.command_line_text
        );
      }
    }
  }

  #[test]
  fn test_config_code_timetree_command_line_and_yaml() {
    let mut config = AppCommand::Timetree.default_config().unwrap();
    config.extend(
      json!({
        "tree": "my data/tree.nwk",
        "metadata": "metadata.tsv",
        "clock_rate": 0.0008,
        "relax": [1.0, 0.5],
        "output_selection": ["auspice", "tracelog"],
        "output_all": "out",
      })
      .as_object()
      .unwrap()
      .clone(),
    );
    let code = config_code(AppCommand::Timetree, &config).unwrap();
    assert_eq!(
      (
        indoc! {r#"
          treetime timetree \
            --tree 'my data/tree.nwk' \
            --metadata metadata.tsv \
            --clock-rate 0.0008 \
            --relax 1.0 0.5 \
            --output-selection 'auspice,tracelog' \
            --output-all out"#}
        .to_owned(),
        indoc! {r#"
          # yaml-language-server: $schema=https://raw.githubusercontent.com/neherlab/treetime/rust/packages/schemas/input-config-timetree.schema.json
          # treetime timetree --config run.yaml
          tree: "my data/tree.nwk"
          metadata: "metadata.tsv"
          clock_rate: 0.0008
          relax:
            - 1.0
            - 0.5
          output_selection:
            - "auspice"
            - "tracelog"
          output_all: "out"
        "#}
        .to_owned(),
      ),
      (code.command_line_text, code.yaml_text)
    );
  }

  #[test]
  fn test_config_code_line_kinds() {
    let mut config = AppCommand::Prune.default_config().unwrap();
    config.extend(
      json!({ "tree": "t.nwk", "prune_empty": true })
        .as_object()
        .unwrap()
        .clone(),
    );
    let code = config_code(AppCommand::Prune, &config).unwrap();
    assert_eq!(
      vec![
        (CodeLineKind::Command, "treetime prune"),
        (CodeLineKind::Input, "--tree t.nwk"),
        (CodeLineKind::Changed, "--prune-empty"),
        (CodeLineKind::Output, "--output-all out"),
      ],
      code
        .command_line
        .iter()
        .map(|line| (line.kind, line.text.as_str()))
        .collect_vec()
    );
  }

  #[test]
  fn test_config_code_empty_list_has_no_command_line_form() {
    let mut config = AppCommand::Clock.default_config().unwrap();
    config.insert("metadata_id_columns".to_owned(), json!([]));
    let code = config_code(AppCommand::Clock, &config).unwrap();
    assert_eq!(
      (
        "# metadata_id_columns = [] has no command-line form; use the YAML config",
        CodeLineKind::Comment,
        "metadata_id_columns: []"
      ),
      (
        code.command_line[1].text.as_str(),
        code.command_line[1].kind,
        code.yaml[2].text.as_str()
      )
    );
  }

  #[test]
  fn test_config_code_nested_setting_is_nested_in_yaml() {
    let mut config = AppCommand::Clock.default_config().unwrap();
    setting_mut(&mut config, &["branch_split".to_owned(), "method".to_owned()])
      .map(|slot| *slot = json!("brent"))
      .unwrap();
    let code = config_code(AppCommand::Clock, &config).unwrap();
    assert_eq!(
      vec!["branch_split:", "  method: \"brent\""],
      code
        .yaml
        .iter()
        .filter(|line| line.kind == CodeLineKind::Changed)
        .map(|line| line.text.as_str())
        .collect_vec()
    );
  }

  #[test]
  fn test_config_code_yaml_of_each_setting_loads_back_to_the_same_config() {
    for command in AppCommand::iter() {
      let defaults = command.default_config().unwrap();
      for spec in command_settings(command).unwrap().settings {
        if spec.role == SettingRole::Output {
          continue;
        }
        let mut config = defaults.clone();
        set_setting(&mut config, &spec.path, changed_value(&spec));
        let code = config_code(command, &config).unwrap();
        assert_eq!(
          Value::Object(with_output_dir(config)),
          load_yaml(command, &code.yaml_text),
          "{command}: `{}`\n{}",
          spec.key,
          code.yaml_text
        );
      }
    }
  }

  #[test]
  fn test_config_code_yaml_passes_the_config_check() {
    let mut config = AppCommand::Ancestral.default_config().unwrap();
    config.extend(
      json!({ "tree": "t.nwk", "alignment": ["a.fasta", "b c.fasta"], "gap_fill": "all", "dense": true })
        .as_object()
        .unwrap()
        .clone(),
    );
    let code = config_code(AppCommand::Ancestral, &config).unwrap();
    let loaded = AppCommand::Ancestral.prepare_text("run.yaml", &code.yaml_text).unwrap();
    assert_eq!(Value::Object(with_output_dir(config)), Value::Object(loaded.config));
  }

  #[test]
  fn test_config_code_output_dir_defaults_to_the_run_folder() {
    let config: Map<String, Value> = AppCommand::Prune.default_config().unwrap();
    let code = config_code(AppCommand::Prune, &config).unwrap();
    assert_eq!("treetime prune \\\n  --output-all out", code.command_line_text);
  }

  mod helpers {
    use crate::command::AppCommand;
    use crate::commands::ancestral::args::TreetimeAncestralArgsRaw;
    use crate::commands::clock::args::TreetimeClockArgsRaw;
    use crate::commands::mugration::args::TreetimeMugrationArgsRaw;
    use crate::commands::optimize::args::TreetimeOptimizeArgsRaw;
    use crate::commands::prune::args::TreetimePruneArgsRaw;
    use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
    use crate::config::catalog::{ListItemKind, SettingKind, SettingRole, SettingSpec};
    use crate::config::load::load_config_document;
    use crate::config::source::ConfigSource;
    use clap::FromArgMatches;
    use serde::Serialize;
    use serde_json::{Map, Value, json};

    fn default_of(spec: &SettingSpec) -> Option<&Value> {
      spec.default_value.as_deref()
    }

    pub(super) fn changed_value(spec: &SettingSpec) -> Value {
      match (spec.role, spec.kind) {
        (SettingRole::Input | SettingRole::InputTemplate, SettingKind::List) => json!(["a b.fasta", "c.fasta"]),
        (SettingRole::Input | SettingRole::InputTemplate, _) => json!("my data/file one.txt"),
        (_, SettingKind::Switch) => json!(true),
        (_, SettingKind::Tristate) => json!(!default_of(spec).and_then(Value::as_bool).unwrap_or(false)),
        (_, SettingKind::Enum) => spec
          .options
          .iter()
          .map(|option| json!(option.value))
          .find(|value| Some(value) != default_of(spec))
          .unwrap(),
        (_, SettingKind::EnumList) => json!(
          spec
            .options
            .iter()
            .map(|option| option.value.clone())
            .filter(|value| value != "All")
            .take(2)
            .collect::<Vec<_>>()
        ),
        (_, SettingKind::Integer) => json!(default_of(spec).and_then(Value::as_u64).unwrap_or(0) + 3),
        (_, SettingKind::Number) => json!(default_of(spec).and_then(Value::as_f64).unwrap_or(0.0) + 0.5),
        (_, SettingKind::Text) => json!("x"),
        (_, SettingKind::List) => match spec.item_kind {
          ListItemKind::String => json!(["x", "y"]),
          ListItemKind::Number => json!([1.5, 2.5]),
          ListItemKind::Integer => json!([1, 2]),
        },
      }
    }

    pub(super) fn set_setting(config: &mut Map<String, Value>, path: &[String], value: Value) {
      let (last, parents) = path.split_last().unwrap();
      let parent = parents.iter().fold(config, |map, key| {
        map
          .entry(key.clone())
          .or_insert_with(|| json!({}))
          .as_object_mut()
          .unwrap()
      });
      parent.insert(last.clone(), value);
    }

    pub(super) fn with_output_dir(mut config: Map<String, Value>) -> Map<String, Value> {
      config.insert("output_all".to_owned(), json!("out"));
      config
    }

    pub(super) fn parse_command_line(command: AppCommand, text: &str) -> Map<String, Value> {
      let words = shlex::split(&text.replace("\\\n", " ")).unwrap();
      let argv = words.into_iter().skip(1).collect::<Vec<_>>();
      let matches = command.cli_command().try_get_matches_from(argv).unwrap();
      match command {
        AppCommand::Timetree => settings(&TreetimeTimetreeArgsRaw::from_arg_matches(&matches).unwrap()),
        AppCommand::Optimize => settings(&TreetimeOptimizeArgsRaw::from_arg_matches(&matches).unwrap()),
        AppCommand::Prune => settings(&TreetimePruneArgsRaw::from_arg_matches(&matches).unwrap()),
        AppCommand::Ancestral => settings(&TreetimeAncestralArgsRaw::from_arg_matches(&matches).unwrap()),
        AppCommand::Clock => settings(&TreetimeClockArgsRaw::from_arg_matches(&matches).unwrap()),
        AppCommand::Mugration => settings(&TreetimeMugrationArgsRaw::from_arg_matches(&matches).unwrap()),
      }
    }

    pub(super) fn load_yaml(command: AppCommand, text: &str) -> Value {
      let source = ConfigSource::new("run.yaml", text);
      match command {
        AppCommand::Timetree => load_config_document::<TreetimeTimetreeArgsRaw>(&source, text, None),
        AppCommand::Optimize => load_config_document::<TreetimeOptimizeArgsRaw>(&source, text, None),
        AppCommand::Prune => load_config_document::<TreetimePruneArgsRaw>(&source, text, None),
        AppCommand::Ancestral => load_config_document::<TreetimeAncestralArgsRaw>(&source, text, None),
        AppCommand::Clock => load_config_document::<TreetimeClockArgsRaw>(&source, text, None),
        AppCommand::Mugration => load_config_document::<TreetimeMugrationArgsRaw>(&source, text, None),
      }
      .unwrap()
    }

    fn settings(raw: &impl Serialize) -> Map<String, Value> {
      serde_json::to_value(raw).unwrap().as_object().unwrap().clone()
    }
  }
}
