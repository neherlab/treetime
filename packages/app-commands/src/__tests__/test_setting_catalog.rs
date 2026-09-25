#[cfg(test)]
mod tests {
  use crate::check_inputs::{InputKind, InputNeed, InputSlot};
  use crate::command::AppCommand;
  use crate::config::catalog::{
    ListItemKind, SettingKind, SettingOption, SettingRole, SettingSpec, command_settings, setting_catalog,
  };
  use helpers::spec;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use strum::IntoEnumIterator;
  use treetime_utils::{o, vec_of_owned};

  #[test]
  fn test_setting_catalog_lists_every_app_command_once() {
    let commands = setting_catalog()
      .unwrap()
      .commands
      .into_iter()
      .map(|settings| settings.command)
      .collect_vec();
    assert_eq!(AppCommand::iter().collect_vec(), commands);
  }

  #[test]
  fn test_setting_catalog_timetree_groups_follow_the_help_headings() {
    assert_eq!(
      vec_of_owned![
        "Input data",
        "Molecular clock",
        "Branch lengths",
        "Dating",
        "Polytomies",
        "Coalescent prior",
        "Plots",
        "Rooting",
        "Substitution model",
        "Ancestral reconstruction",
        "Output",
        "Tree ordering",
        "Reproducibility",
      ],
      command_settings(AppCommand::Timetree).unwrap().groups
    );
  }

  #[test]
  fn test_setting_catalog_settings_are_ordered_by_group() {
    for command in AppCommand::iter() {
      let settings = command_settings(command).unwrap();
      let group_order = settings
        .settings
        .iter()
        .map(|spec| settings.groups.iter().position(|group| *group == spec.group).unwrap())
        .collect_vec();
      assert!(group_order.is_sorted(), "{command}: {group_order:?}");
    }
  }

  #[test]
  fn test_setting_catalog_has_every_leaf_of_the_config_once() {
    for command in AppCommand::iter() {
      let keys = command_settings(command)
        .unwrap()
        .settings
        .into_iter()
        .map(|spec| spec.key)
        .sorted()
        .collect_vec();
      let unique = keys.iter().unique().count();
      assert_eq!(keys.len(), unique, "{command}");
    }
  }

  #[test]
  fn test_setting_catalog_clock_rate() {
    assert_eq!(
      SettingSpec {
        key: o!("clock_rate"),
        path: vec_of_owned!["clock_rate"],
        label: o!("Clock rate"),
        flag: o!("--clock-rate"),
        group: o!("Molecular clock"),
        role: SettingRole::Setting,
        kind: SettingKind::Number,
        nullable: true,
        options: vec![],
        item_kind: ListItemKind::String,
        default_value: Value::Null,
        minimum: None,
        examples: vec![json!(0.001)],
        value_names: vec_of_owned!["CLOCK_RATE"],
        conflicts: vec![],
        help: o!("If specified, the rate of the molecular clock won't be optimized."),
        more: o!(""),
      },
      spec(AppCommand::Timetree, "clock_rate")
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::declared_on_the_setting( AppCommand::Timetree, "coalescent_opt",     vec_of_owned!["coalescent", "coalescent_skyline"])]
  #[case::declared_on_the_others(  AppCommand::Timetree, "coalescent",         vec_of_owned!["coalescent_opt", "coalescent_skyline"])]
  #[case::keep_root(               AppCommand::Clock,    "keep_root",          vec_of_owned!["reroot", "reroot_tips"])]
  #[case::reroot(                  AppCommand::Clock,    "reroot",             vec_of_owned!["keep_root", "reroot_tips"])]
  #[case::none(                    AppCommand::Timetree, "max_iter",           vec![])]
  #[trace]
  fn test_setting_catalog_conflicts_are_symmetric(
    #[case] command: AppCommand,
    #[case] key: &str,
    #[case] expected: Vec<String>,
  ) {
    assert_eq!(expected, spec(command, key).conflicts);
  }

  #[test]
  fn test_setting_catalog_relax_names_its_values_and_suggests_a_weak_prior() {
    let relax = spec(AppCommand::Timetree, "relax");
    assert_eq!(
      (vec_of_owned!["SLACK", "COUPLING"], vec![json!([1.0, 0.0])]),
      (relax.value_names, relax.examples)
    );
  }

  #[test]
  fn test_setting_catalog_description_splits_into_help_and_more() {
    let time_marginal = spec(AppCommand::Timetree, "time_marginal");
    assert_eq!(
      (
        "Control when marginal time distributions are used for output.",
        "All modes use marginal inference during optimization. The mode controls whether confidence intervals are extracted from the resulting distributions:\n\n- `never`: no confidence interval output (default) - `always`: write confidence intervals from distributions computed during optimization - `only-final`: run one extra inference pass after optimization, then write confidence intervals"
      ),
      (time_marginal.help.as_str(), time_marginal.more.as_str())
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::switch(         AppCommand::Timetree,  "confidence",          SettingKind::Switch,   false, json!(false))]
  #[case::tristate(       AppCommand::Ancestral, "dense",               SettingKind::Tristate, true,  Value::Null)]
  #[case::integer(        AppCommand::Timetree,  "max_iter",            SettingKind::Integer,  false, json!(2))]
  #[case::optional_int(   AppCommand::Timetree,  "seed",                SettingKind::Integer,  true,  Value::Null)]
  #[case::text(           AppCommand::Mugration, "missing_data",        SettingKind::Text,     false, json!("?"))]
  #[case::number_list(    AppCommand::Timetree,  "relax",               SettingKind::List,     false, json!([]))]
  #[case::string_list(    AppCommand::Clock,     "metadata_id_columns", SettingKind::List,     false, json!(["strain", "name", "accession"]))]
  #[case::enum_list(      AppCommand::Clock,     "output_selection",    SettingKind::EnumList, false, json!([]))]
  #[case::nested_enum(    AppCommand::Clock,     "branch_split.method", SettingKind::Enum,     false, json!("grid"))]
  #[case::nested_integer( AppCommand::Clock,     "branch_split.n_points", SettingKind::Integer, false, json!(11))]
  #[trace]
  fn test_setting_catalog_kind_and_default(
    #[case] command: AppCommand,
    #[case] key: &str,
    #[case] kind: SettingKind,
    #[case] nullable: bool,
    #[case] default_value: Value,
  ) {
    let spec = spec(command, key);
    assert_eq!((kind, nullable, default_value), (spec.kind, spec.nullable, spec.default_value));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::setting(       AppCommand::Timetree,  "clock_filter",       SettingRole::Setting)]
  #[case::input(         AppCommand::Timetree,  "tree",               SettingRole::Input)]
  #[case::input_list(    AppCommand::Ancestral, "alignment",          SettingRole::Input)]
  #[case::input_template(AppCommand::Ancestral, "translations",       SettingRole::InputTemplate)]
  #[case::output(        AppCommand::Timetree,  "output_clock_model", SettingRole::Output)]
  #[case::plot(          AppCommand::Clock,     "plot_rtt",           SettingRole::Output)]
  #[trace]
  fn test_setting_catalog_role(#[case] command: AppCommand, #[case] key: &str, #[case] role: SettingRole) {
    assert_eq!(role, spec(command, key).role);
  }

  #[test]
  fn test_setting_catalog_list_items_of_relax_are_numbers() {
    assert_eq!(ListItemKind::Number, spec(AppCommand::Timetree, "relax").item_kind);
  }

  #[test]
  fn test_setting_catalog_optional_enum_lists_its_values() {
    let reroot = spec(AppCommand::Timetree, "reroot");
    assert_eq!(
      (
        SettingKind::Enum,
        vec_of_owned!["least-squares", "min-dev", "oldest", "clock-filter"]
      ),
      (
        reroot.kind,
        reroot.options.into_iter().map(|option| option.value).collect_vec()
      )
    );
  }

  #[test]
  fn test_setting_catalog_enum_options_carry_their_descriptions() {
    let options = spec(AppCommand::Timetree, "time_marginal").options;
    assert_eq!(
      vec![
        SettingOption {
          value: o!("never"),
          help: o!("")
        },
        SettingOption {
          value: o!("always"),
          help: o!("")
        },
        SettingOption {
          value: o!("only-final"),
          help: o!("")
        },
      ],
      options
    );
  }

  #[test]
  fn test_setting_catalog_every_enum_value_has_a_command_line_spelling() {
    for command in AppCommand::iter() {
      let cli = command.cli_command();
      for spec in command_settings(command).unwrap().settings {
        let arg = cli
          .get_arguments()
          .find(|arg| arg.get_long() == spec.flag.strip_prefix("--"))
          .unwrap();
        let spellings = arg
          .get_possible_values()
          .iter()
          .map(|value| value.get_name().replace(['-', '_'], "").to_lowercase())
          .collect_vec();
        for option in &spec.options {
          assert!(
            spellings.contains(&option.value.replace(['-', '_'], "").to_lowercase()),
            "{command}: `{}` value `{}` has no command-line spelling among {spellings:?}",
            spec.key,
            option.value
          );
        }
      }
    }
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree( AppCommand::Timetree,  true,  &[(InputKind::Tree, InputNeed::Required), (InputKind::Metadata, InputNeed::Required), (InputKind::Alignment, InputNeed::Recommended)])]
  #[case::clock(    AppCommand::Clock,     true,  &[(InputKind::Tree, InputNeed::Required), (InputKind::Metadata, InputNeed::Required), (InputKind::Alignment, InputNeed::Optional)])]
  #[case::ancestral(AppCommand::Ancestral, false, &[(InputKind::Tree, InputNeed::Required), (InputKind::Alignment, InputNeed::Required)])]
  #[case::mugration(AppCommand::Mugration, false, &[(InputKind::Tree, InputNeed::Required), (InputKind::Metadata, InputNeed::Required)])]
  #[case::optimize( AppCommand::Optimize,  false, &[(InputKind::Tree, InputNeed::Required), (InputKind::Alignment, InputNeed::Required)])]
  #[case::prune(    AppCommand::Prune,     false, &[(InputKind::Tree, InputNeed::Required), (InputKind::Alignment, InputNeed::Optional)])]
  #[trace]
  fn test_setting_catalog_inputs_of_command(
    #[case] command: AppCommand,
    #[case] uses_dates: bool,
    #[case] inputs: &[(InputKind, InputNeed)],
  ) {
    let settings = command_settings(command).unwrap();
    let expected = inputs.iter().map(|(kind, need)| (*kind, *need)).collect_vec();
    let actual = settings.inputs.iter().map(|input| (input.kind, input.need)).collect_vec();
    assert_eq!((uses_dates, expected), (settings.uses_dates, actual));
  }

  #[test]
  fn test_setting_catalog_input_slots_describe_the_readers() {
    let settings = command_settings(AppCommand::Timetree).unwrap();
    assert_eq!(
      vec![
        InputSlot {
          kind: InputKind::Tree,
          need: InputNeed::Required,
          label: o!("Tree"),
          formats: o!("Newick"),
          extensions: vec_of_owned!["nwk", "newick", "tree", "tre", "bz2", "xz", "zst", "gz"],
          list: false,
        },
        InputSlot {
          kind: InputKind::Metadata,
          need: InputNeed::Required,
          label: o!("Metadata"),
          formats: o!("CSV, TSV or SSV table"),
          extensions: vec_of_owned!["csv", "tsv", "ssv", "bz2", "xz", "zst", "gz"],
          list: false,
        },
        InputSlot {
          kind: InputKind::Alignment,
          need: InputNeed::Recommended,
          label: o!("Alignment"),
          formats: o!("Aligned FASTA"),
          extensions: vec_of_owned!["fasta", "fa", "fas", "aln", "bz2", "xz", "zst", "gz"],
          list: true,
        },
      ],
      settings.inputs
    );
  }

  #[test]
  fn test_setting_catalog_every_input_is_a_setting_of_its_command() {
    for command in AppCommand::iter() {
      let settings = command_settings(command).unwrap();
      for input in &settings.inputs {
        let spec = settings
          .settings
          .iter()
          .find(|spec| spec.key == input.kind.setting())
          .unwrap_or_else(|| panic!("{command}: no setting `{}`", input.kind.setting()));
        assert_eq!(SettingRole::Input, spec.role, "{command}: {}", spec.key);
      }
    }
  }

  mod helpers {
    use crate::command::AppCommand;
    use crate::config::catalog::{SettingSpec, command_settings};

    pub(super) fn spec(command: AppCommand, key: &str) -> SettingSpec {
      command_settings(command)
        .unwrap()
        .settings
        .into_iter()
        .find(|spec| spec.key == key)
        .unwrap()
    }
  }
}
