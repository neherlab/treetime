#[cfg(test)]
mod tests {
  use crate::check_inputs::{InputFacts, InputKind, InputProblem};
  use crate::command::AppCommand;
  use crate::json_value::JsonValue;
  use crate::run_checks::{
    CheckContext, CheckFix, CheckLevel, ConfigRejection, RunCheck, SettingPatch, rejection_messages, run_checks,
  };
  use helpers::{checks, config, dates, facts_with_tree, metadata, problem, texts};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use treetime_utils::{o, vec_of_owned};

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree_needs_tree_and_metadata(AppCommand::Timetree,  json!({}),                  &["missing-tree", "missing-metadata"])]
  #[case::timetree_alignment_recommended( AppCommand::Timetree,  json!({ "tree": "t.nwk", "metadata": "m.tsv" }), &[])]
  #[case::ancestral_needs_alignment(      AppCommand::Ancestral, json!({ "tree": "t.nwk" }), &["missing-alignment"])]
  #[case::empty_alignment_list_is_missing(AppCommand::Optimize,  json!({ "tree": "t.nwk", "alignment": [] }), &["missing-alignment"])]
  #[case::empty_path_is_missing(          AppCommand::Prune,     json!({ "tree": "" }),      &["missing-tree"])]
  #[case::prune_alignment_optional(       AppCommand::Prune,     json!({ "tree": "t.nwk" }), &[])]
  #[trace]
  fn test_run_checks_missing_inputs(#[case] command: AppCommand, #[case] settings: Value, #[case] ids: &[&str]) {
    let config = config(command, &settings);
    let found = checks(command, Some(&config), None, None);
    assert_eq!(
      ids.iter().map(|id| ((*id).to_owned(), CheckLevel::Block)).collect::<Vec<_>>(),
      found.into_iter().map(|check| (check.id, check.level)).collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_run_checks_missing_input_texts() {
    let config = config(AppCommand::Ancestral, &json!({}));
    assert_eq!(
      vec!["Add a tree file.", "Add an alignment file."],
      texts(&checks(AppCommand::Ancestral, Some(&config), None, None))
    );
  }

  #[test]
  fn test_run_checks_unreadable_input_blocks() {
    let config = config(AppCommand::Prune, &json!({ "tree": "t.nwk" }));
    let facts = InputFacts {
      problems: vec![InputProblem {
        input: InputKind::Tree,
        message: o!("When reading 't.nwk': no such file"),
      }],
      ..InputFacts::default()
    };
    assert_eq!(
      vec![RunCheck {
        id: o!("unreadable-tree"),
        level: CheckLevel::Block,
        text: o!("The tree cannot be read: When reading 't.nwk': no such file"),
        settings: vec_of_owned!["tree"],
        fix: None,
      }],
      checks(AppCommand::Prune, Some(&config), None, Some(&facts))
    );
  }

  #[test]
  fn test_run_checks_rejected_config_lists_each_problem_with_its_help() {
    let config = config(AppCommand::Prune, &json!({ "tree": "t.nwk" }));
    let problems = vec![
      problem("unknown field `prune_shrot`", Some("did you mean `prune_short`?")),
      problem("invalid type: string \"x\", expected a number", None),
    ];
    let rejection = ConfigRejection {
      message: "invalid configuration: ...",
      causes: &[],
      problems: &problems,
    };
    assert_eq!(
      vec![
        (
          o!("config-0"),
          o!("unknown field `prune_shrot` (did you mean `prune_short`?)")
        ),
        (o!("config-1"), o!("invalid type: string \"x\", expected a number")),
      ],
      checks(AppCommand::Prune, Some(&config), Some(rejection), None)
        .into_iter()
        .map(|check| (check.id, check.text))
        .collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_run_checks_rejected_config_without_problems_shows_the_error_and_its_causes() {
    let config = config(AppCommand::Prune, &json!({ "tree": "t.nwk" }));
    let causes = vec_of_owned!["the value is negative"];
    let rejection = ConfigRejection {
      message: "invalid --prune-short",
      causes: &causes,
      problems: &[],
    };
    assert_eq!(
      vec!["invalid --prune-short: the value is negative"],
      texts(&checks(AppCommand::Prune, Some(&config), Some(rejection), None))
    );
  }

  #[test]
  fn test_run_checks_rejection_that_repeats_a_missing_input_is_left_out() {
    let config = config(AppCommand::Prune, &json!({}));
    let rejection = ConfigRejection {
      message: "the following required arguments were not provided:\n  --tree <TREE>",
      causes: &[],
      problems: &[],
    };
    assert_eq!(
      vec!["Add a tree file."],
      texts(&checks(AppCommand::Prune, Some(&config), Some(rejection), None))
    );
  }

  #[test]
  fn test_run_checks_rejection_messages_without_problems_only() {
    let causes = vec_of_owned!["b", "c"];
    let rejection = ConfigRejection {
      message: "a",
      causes: &causes,
      problems: &[],
    };
    assert_eq!(
      (vec_of_owned!["a: b: c"], Vec::<String>::new()),
      (
        rejection_messages(&rejection, false),
        rejection_messages(&rejection, true)
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::reads_alignment(  AppCommand::Ancestral, vec!["tips-without-sequence: 5 of 20 tree tips have no sequence in the alignment: A, B, C, and 2 more"])]
  #[case::without_alignment(AppCommand::Mugration, vec![])]
  #[trace]
  fn test_run_checks_tips_without_sequence(#[case] command: AppCommand, #[case] expected: Vec<&str>) {
    let config = config(command, &json!({ "tree": "t.nwk", "alignment": ["a.fasta"], "metadata": "m.tsv" }));
    let facts = InputFacts {
      tips_without_sequence: Some(vec_of_owned!["A", "B", "C", "D", "E"]),
      ..facts_with_tree(20)
    };
    assert_eq!(
      expected,
      checks(command, Some(&config), None, Some(&facts))
        .iter()
        .map(|check| format!("{}: {}", check.id, check.text))
        .collect::<Vec<_>>()
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dates(       AppCommand::Timetree,  "2 of 20 tree tips have no metadata row and get no date: A, B")]
  #[case::traits(      AppCommand::Mugration, "2 of 20 tree tips have no metadata row and get no trait value: A, B")]
  #[trace]
  fn test_run_checks_tips_without_metadata_warn(#[case] command: AppCommand, #[case] expected: &str) {
    let config = config(command, &json!({ "tree": "t.nwk", "metadata": "m.tsv" }));
    let facts = InputFacts {
      tips_without_metadata: Some(vec_of_owned!["A", "B"]),
      ..facts_with_tree(20)
    };
    let found = checks(command, Some(&config), None, Some(&facts));
    assert_eq!(
      vec![(o!("tips-without-metadata"), CheckLevel::Warn, expected.to_owned())],
      found.into_iter().map(|check| (check.id, check.level, check.text)).collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_run_checks_tips_without_metadata_ignored_by_commands_without_metadata() {
    let config = config(
      AppCommand::Optimize,
      &json!({ "tree": "t.nwk", "alignment": ["a.fasta"] }),
    );
    let facts = InputFacts {
      tips_without_metadata: Some(vec_of_owned!["A"]),
      ..facts_with_tree(20)
    };
    assert_eq!(
      Vec::<RunCheck>::new(),
      checks(AppCommand::Optimize, Some(&config), None, Some(&facts))
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree( AppCommand::Timetree,  vec!["no-date-column: The metadata has no date column. Set the date column to one of: strain, when, country."])]
  #[case::clock(    AppCommand::Clock,     vec!["no-date-column: The metadata has no date column. Set the date column to one of: strain, when, country."])]
  #[case::mugration(AppCommand::Mugration, vec![])]
  #[trace]
  fn test_run_checks_no_date_column(#[case] command: AppCommand, #[case] expected: Vec<&str>) {
    let config = config(command, &json!({ "tree": "t.nwk", "metadata": "m.tsv" }));
    let facts = InputFacts {
      metadata: Some(metadata(None, None)),
      ..facts_with_tree(20)
    };
    assert_eq!(
      expected,
      checks(command, Some(&config), None, Some(&facts))
        .iter()
        .map(|check| format!("{}: {}", check.id, check.text))
        .collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_run_checks_unreadable_dates_warn() {
    let config = config(AppCommand::Clock, &json!({ "tree": "t.nwk", "metadata": "m.tsv" }));
    let facts = InputFacts {
      metadata: Some(metadata(Some("date"), Some(dates(vec_of_owned!["A", "B"], 0, 0)))),
      ..facts_with_tree(20)
    };
    assert_eq!(
      vec![RunCheck {
        id: o!("unreadable-dates"),
        level: CheckLevel::Warn,
        text: o!("2 dates cannot be read and those samples get no date: A, B. Use 2015-06-21, 2015-06-XX or 2015.47."),
        settings: vec![],
        fix: None,
      }],
      checks(AppCommand::Clock, Some(&config), None, Some(&facts))
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::without_rate_uncertainty(AppCommand::Timetree, json!({ "confidence": true }),                              true)]
  #[case::with_covariation(        AppCommand::Timetree, json!({ "confidence": true, "covariation": true }),         false)]
  #[case::with_clock_std_dev(      AppCommand::Timetree, json!({ "confidence": true, "clock_std_dev": 0.0001 }),     false)]
  #[case::without_confidence(      AppCommand::Timetree, json!({}),                                                  false)]
  #[trace]
  fn test_run_checks_confidence_without_rate_uncertainty(
    #[case] command: AppCommand,
    #[case] settings: Value,
    #[case] warns: bool,
  ) {
    let mut settings = settings;
    settings["tree"] = json!("t.nwk");
    settings["metadata"] = json!("m.tsv");
    let config = config(command, &settings);
    let expected = warns.then(|| RunCheck {
      id: o!("confidence-without-rate-uncertainty"),
      level: CheckLevel::Warn,
      text: o!("Date intervals need rate uncertainty: without the covariation-aware regression or a clock rate std. dev., this run writes no intervals."),
      settings: vec_of_owned!["confidence", "covariation", "clock_std_dev"],
      fix: Some(CheckFix {
        label: o!("Use covariation"),
        patch: vec![SettingPatch { path: vec_of_owned!["covariation"], value: Some(JsonValue(json!(true))) }],
      }),
    });
    assert_eq!(expected.into_iter().collect::<Vec<_>>(), checks(command, Some(&config), None, None));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::above_a_quarter(  AppCommand::Timetree,  3, 10, true)]
  #[case::exactly_a_quarter(AppCommand::Timetree,  2,  8, false)]
  #[case::no_exact_days(    AppCommand::Clock,     0,  0, false)]
  #[case::without_dates(    AppCommand::Mugration, 5,  5, false)]
  #[trace]
  fn test_run_checks_dates_rounded_to_the_month(
    #[case] command: AppCommand,
    #[case] on_day_1_or_15: usize,
    #[case] exact_days: usize,
    #[case] advises: bool,
  ) {
    let config = config(command, &json!({ "tree": "t.nwk", "metadata": "m.tsv" }));
    let facts = InputFacts {
      metadata: Some(metadata(Some("date"), Some(dates(vec![], exact_days, on_day_1_or_15)))),
      ..facts_with_tree(20)
    };
    let expected = advises.then(|| RunCheck {
      id: o!("dates-on-day-1-or-15"),
      level: CheckLevel::Advice,
      text: format!("{on_day_1_or_15} of {exact_days} dates fall on the 1st or 15th of a month. If only the month is known, write it as 2015-06-XX so TreeTime treats it as a range."),
      settings: vec![],
      fix: None,
    });
    assert_eq!(expected.into_iter().collect::<Vec<_>>(), checks(command, Some(&config), None, Some(&facts)));
  }

  #[test]
  fn test_run_checks_are_listed_in_rule_order() {
    let config = config(
      AppCommand::Timetree,
      &json!({ "metadata": "m.tsv", "confidence": true }),
    );
    let facts = InputFacts {
      tips_without_metadata: Some(vec_of_owned!["A"]),
      metadata: Some(metadata(Some("date"), Some(dates(vec_of_owned!["B"], 4, 4)))),
      problems: vec![InputProblem {
        input: InputKind::Alignment,
        message: o!("bad"),
      }],
      ..facts_with_tree(20)
    };
    let problems = vec![problem("unknown field `x`", None)];
    let rejection = ConfigRejection {
      message: "invalid configuration",
      causes: &[],
      problems: &problems,
    };
    assert_eq!(
      vec_of_owned![
        "missing-tree",
        "unreadable-alignment",
        "config-0",
        "tips-without-metadata",
        "unreadable-dates",
        "confidence-without-rate-uncertainty",
        "dates-on-day-1-or-15"
      ],
      checks(AppCommand::Timetree, Some(&config), Some(rejection), Some(&facts))
        .into_iter()
        .map(|check| check.id)
        .collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_run_checks_without_parsed_settings_check_only_the_rejection() {
    let rejection = ConfigRejection {
      message: "invalid configuration: could not parse config",
      causes: &[],
      problems: &[],
    };
    assert_eq!(
      vec!["invalid configuration: could not parse config"],
      texts(&run_checks(&CheckContext {
        command: AppCommand::Prune,
        config: None,
        rejection: Some(rejection),
        facts: None,
      }))
    );
  }

  mod helpers {
    use crate::check_inputs::{DateFacts, MetadataFacts};
    use crate::check_inputs::{InputFacts, TreeFacts};
    use crate::command::AppCommand;
    use crate::config::source::ConfigProblem;
    use crate::run_checks::{CheckContext, ConfigRejection, RunCheck, run_checks};
    use serde_json::{Map, Value};
    use treetime_utils::{o, vec_of_owned};

    pub(super) fn config(command: AppCommand, settings: &Value) -> Map<String, Value> {
      command.config_over_defaults(settings).unwrap()
    }

    pub(super) fn checks(
      command: AppCommand,
      config: Option<&Map<String, Value>>,
      rejection: Option<ConfigRejection<'_>>,
      facts: Option<&InputFacts>,
    ) -> Vec<RunCheck> {
      run_checks(&CheckContext {
        command,
        config,
        rejection,
        facts,
      })
    }

    pub(super) fn texts(checks: &[RunCheck]) -> Vec<&str> {
      checks.iter().map(|check| check.text.as_str()).collect()
    }

    pub(super) fn facts_with_tree(tips: usize) -> InputFacts {
      InputFacts {
        tree: Some(TreeFacts {
          tips,
          internal_nodes: tips - 1,
          polytomies: 0,
          unnamed_tips: 0,
          duplicate_tip_names: vec![],
        }),
        ..InputFacts::default()
      }
    }
    pub(super) fn problem(message: &str, help: Option<&str>) -> ConfigProblem {
      ConfigProblem {
        code: o!("config::test"),
        message: message.to_owned(),
        span: None,
        help: help.map(str::to_owned),
      }
    }

    pub(super) fn metadata(date_column: Option<&str>, dates: Option<DateFacts>) -> MetadataFacts {
      MetadataFacts {
        rows: 20,
        columns: vec_of_owned!["strain", "when", "country"],
        id_column: o!("strain"),
        date_column: date_column.map(str::to_owned),
        dates,
      }
    }

    pub(super) fn dates(unreadable: Vec<String>, exact_days: usize, on_day_1_or_15: usize) -> DateFacts {
      DateFacts {
        readable: 20 - unreadable.len(),
        unreadable,
        exact_days,
        on_day_1_or_15,
      }
    }
  }
}
