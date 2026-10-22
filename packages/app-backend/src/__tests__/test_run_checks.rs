#[cfg(test)]
mod tests {
  use crate::check_inputs::{AlignmentFacts, InputFacts, InputKind, InputProblem};
  use crate::run_checks::{
    CheckContext, CheckFix, CheckLevel, ConfigRejection, RunCheck, SettingPatch, rejection_messages, run_checks,
  };
  use app_commands::command::AppCommand;
  use app_commands::json_value::JsonValue;
  use helpers::{alignment, checks, config, dates, facts_with_tree, metadata, problem, texts};
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
        partition: None,
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
      problem("unknown field `collapse_shrot`", Some("did you mean `collapse_short`?")),
      problem("invalid type: string \"x\", expected a number", None),
    ];
    let rejection = ConfigRejection {
      message: "invalid configuration: ...",
      causes: &[],
      problems: &problems,
      settings: &[],
    };
    assert_eq!(
      vec![
        (
          o!("config-0"),
          o!("unknown field `collapse_shrot` (did you mean `collapse_short`?)")
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
      message: "invalid --collapse-short",
      causes: &causes,
      problems: &[],
      settings: &[],
    };
    assert_eq!(
      vec!["invalid --collapse-short: the value is negative"],
      texts(&checks(AppCommand::Prune, Some(&config), Some(rejection), None))
    );
  }

  #[test]
  fn test_run_checks_rejection_that_repeats_a_missing_input_is_left_out() {
    let config = config(AppCommand::Prune, &json!({}));
    let rejection = ConfigRejection {
      message: "the following required arguments were not provided:\n  --tree <FILE>",
      causes: &[],
      problems: &[],
      settings: &[],
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
      settings: &[],
    };
    assert_eq!(
      (vec_of_owned!["a: b: c"], Vec::<String>::new()),
      (
        rejection_messages(&rejection, false),
        rejection_messages(&rejection, true)
      )
    );
  }

  #[test]
  fn test_run_checks_duplicate_names_warn_once_per_input() {
    let config = config(
      AppCommand::Ancestral,
      &json!({ "tree": "t.nwk", "alignment": ["a.fasta"] }),
    );
    let mut facts = facts_with_tree(20);
    if let Some(tree) = facts.tree.as_mut() {
      tree.duplicate_node_names = vec_of_owned!["A", "B", "C", "D", "E"];
    }
    facts.alignments = vec![AlignmentFacts {
      duplicate_names: vec_of_owned!["A"],
      ..alignment(1, "a.fasta")
    }];

    let actual = checks(AppCommand::Ancestral, Some(&config), None, Some(&facts));

    assert_eq!(
      vec![
        (
          o!("duplicate-node-names"),
          CheckLevel::Warn,
          o!(
            "The tree gives the same name to more than one node: A, B, C, and 2 more. TreeTime keeps the name of the first such node and writes the others as <name>_2, <name>_3 and so on. Nodes with the same input name still receive the same data from the other inputs."
          ),
          vec_of_owned!["tree"],
        ),
        (
          o!("duplicate-sequence-names"),
          CheckLevel::Warn,
          o!("The alignment has more than one sequence named A. TreeTime uses the first sequence of each name."),
          vec_of_owned!["alignment"],
        ),
      ],
      actual
        .into_iter()
        .map(|check| (check.id, check.level, check.text, check.settings))
        .collect::<Vec<_>>()
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::ancestral_up_to_a_third(  AppCommand::Ancestral, json!({}),                              5,  vec![("tips-without-sequence", CheckLevel::Warn, "5 of 20 tree tips have no sequence in the alignment: A, B, C, and 2 more. TreeTime treats their sequences as unknown.")])]
  #[case::ancestral_over_a_third(   AppCommand::Ancestral, json!({}),                              7,  vec![("tips-without-sequence", CheckLevel::Block, "7 of 20 tree tips have no sequence in the alignment: A, B, C, and 4 more. More than a third of the tips are missing, so the alignment may not belong to the tree.")])]
  #[case::ancestral_over_but_allowed(AppCommand::Ancestral, json!({ "allow_missing_sequences": true }), 7,  vec![("tips-without-sequence", CheckLevel::Warn, "7 of 20 tree tips have no sequence in the alignment: A, B, C, and 4 more. TreeTime treats their sequences as unknown.")])]
  #[case::optimize_over_a_third(    AppCommand::Optimize,  json!({}),                              7,  vec![("tips-without-sequence", CheckLevel::Block, "7 of 20 tree tips have no sequence in the alignment: A, B, C, and 4 more. More than a third of the tips are missing, so the alignment may not belong to the tree.")])]
  #[case::optimize_over_but_allowed(AppCommand::Optimize, json!({ "allow_missing_sequences": true }), 7,  vec![("tips-without-sequence", CheckLevel::Warn, "7 of 20 tree tips have no sequence in the alignment: A, B, C, and 4 more. TreeTime treats their sequences as unknown.")])]
  #[case::prune_up_to_a_third(      AppCommand::Prune,     json!({}),                              1,  vec![("tips-without-sequence", CheckLevel::Warn, "1 of 20 tree tips have no sequence in the alignment: A. TreeTime treats their sequences as unknown.")])]
  #[case::clock_ignores_sequences(  AppCommand::Clock,     json!({}),                              20, vec![])]
  #[case::mugration_without_alignment(AppCommand::Mugration, json!({}),                            5,  vec![])]
  #[trace]
  fn test_run_checks_tips_without_sequence(
    #[case] command: AppCommand,
    #[case] extra: Value,
    #[case] n_missing: usize,
    #[case] expected: Vec<(&str, CheckLevel, &str)>,
  ) {
    let mut settings = json!({ "tree": "t.nwk", "alignment": ["a.fasta"], "metadata": "m.tsv" });
    settings.as_object_mut().unwrap().extend(extra.as_object().unwrap().clone());
    let config = config(command, &settings);
    let names = ('A'..='T').take(n_missing).map(String::from).collect();
    let facts = InputFacts {
      alignments: vec![AlignmentFacts {
        tips_without_sequence: Some(names),
        ..alignment(1, "a.fasta")
      }],
      ..facts_with_tree(20)
    };
    assert_eq!(
      expected
        .into_iter()
        .map(|(id, level, text)| (id.to_owned(), level, text.to_owned()))
        .collect::<Vec<_>>(),
      checks(command, Some(&config), None, Some(&facts))
        .into_iter()
        .filter(|check| check.id.starts_with("tips-without-sequence"))
        .map(|check| (check.id, check.level, check.text))
        .collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_run_checks_ancestral_offers_to_treat_missing_sequences_as_unknown() {
    let config = config(
      AppCommand::Ancestral,
      &json!({ "tree": "t.nwk", "alignment": ["a.fasta"] }),
    );
    let facts = InputFacts {
      alignments: vec![AlignmentFacts {
        tips_without_sequence: Some(vec_of_owned!["A", "B"]),
        ..alignment(1, "a.fasta")
      }],
      ..facts_with_tree(3)
    };
    let found = checks(AppCommand::Ancestral, Some(&config), None, Some(&facts));
    assert_eq!(
      vec![(
        vec_of_owned!["alignment", "allow_missing_sequences"],
        Some(CheckFix {
          label: o!("Treat missing sequences as unknown"),
          patch: vec![SettingPatch {
            key: o!("allow_missing_sequences"),
            value: Some(JsonValue(Value::Bool(true))),
          }],
        })
      )],
      found
        .into_iter()
        .map(|check| (check.settings, check.fix))
        .collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_run_checks_apply_the_alignment_rules_to_each_file() {
    let config = config(
      AppCommand::Optimize,
      &json!({ "tree": "t.nwk", "alignment": ["ha/aligned.fasta", "na/aligned.fasta", "mp/aligned.fasta"] }),
    );
    let facts = InputFacts {
      alignments: vec![
        AlignmentFacts {
          duplicate_names: vec_of_owned!["A"],
          tips_without_sequence: Some(vec![]),
          ..alignment(1, "ha/aligned.fasta")
        },
        AlignmentFacts {
          min_length: 9,
          tips_without_sequence: Some(vec_of_owned!["A", "B"]),
          ..alignment(2, "na/aligned.fasta")
        },
      ],
      problems: vec![InputProblem {
        input: InputKind::Alignment,
        partition: Some(3),
        message: o!("bad"),
      }],
      ..facts_with_tree(3)
    };
    assert_eq!(
      vec![
        (
          o!("unreadable-alignment-3"),
          CheckLevel::Block,
          o!("The alignment cannot be read: bad")
        ),
        (
          o!("duplicate-sequence-names-1"),
          CheckLevel::Warn,
          o!(
            "The alignment 1 'ha/aligned.fasta' has more than one sequence named A. TreeTime uses the first sequence of each name."
          )
        ),
        (
          o!("unaligned-sequences-2"),
          CheckLevel::Block,
          o!("The alignment 2 'na/aligned.fasta' is not aligned: its sequences have lengths 9 to 10.")
        ),
        (
          o!("tips-without-sequence-2"),
          CheckLevel::Block,
          o!(
            "2 of 3 tree tips have no sequence in the alignment 2 'na/aligned.fasta': A, B. More than a third of the tips are missing, so the alignment may not belong to the tree."
          )
        ),
      ],
      checks(AppCommand::Optimize, Some(&config), None, Some(&facts))
        .into_iter()
        .map(|check| (check.id, check.level, check.text))
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
        partition: Some(1),
        message: o!("bad"),
      }],
      ..facts_with_tree(20)
    };
    let problems = vec![problem("unknown field `x`", None)];
    let rejection = ConfigRejection {
      message: "invalid configuration",
      causes: &[],
      problems: &problems,
      settings: &[],
    };
    assert_eq!(
      vec_of_owned![
        "missing-tree",
        "unreadable-alignment",
        "config-0",
        "tips-without-metadata",
        "unreadable-dates",
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
      settings: &[],
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

  #[test]
  fn test_run_checks_partitions_block_entries_that_leave_sites_out() {
    let config = config(
      AppCommand::Ancestral,
      &json!({
        "tree": "t.nwk",
        "alignment": ["a.fasta", "b.fasta"],
        "partition": [{ "alignment": 1, "sites": ["1-4"] }, { "alignment": 2 }],
      }),
    );
    let facts = InputFacts {
      alignments: vec![alignment(1, "a.fasta"), alignment(2, "b.fasta")],
      ..facts_with_tree(3)
    };
    assert_eq!(
      vec![(
        o!("partitions"),
        CheckLevel::Block,
        o!(
          "No partition takes 6 sites of alignment 1, the first is site 5: list them in an entry, or add an entry without `sites`, which takes the remaining sites"
        ),
        vec_of_owned!["partition", "alignment"],
      )],
      checks(AppCommand::Ancestral, Some(&config), None, Some(&facts))
        .into_iter()
        .map(|check| (check.id, check.level, check.text, check.settings))
        .collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_run_checks_partitions_block_a_model_without_its_alphabet() {
    let config = config(
      AppCommand::Optimize,
      &json!({ "tree": "t.nwk", "alignment": ["a.fasta"], "alphabet": "aa", "model": "hky85" }),
    );
    let facts = InputFacts {
      alignments: vec![alignment(1, "a.fasta")],
      ..facts_with_tree(3)
    };
    assert_eq!(
      vec![(
        o!(
          "Partition 1 has the model hky85, which needs the alphabet nuc because it sets rates between nucleotides, but the alphabet is aa"
        ),
        vec_of_owned!["model", "model_params", "alphabet"],
      )],
      checks(AppCommand::Optimize, Some(&config), None, Some(&facts))
        .into_iter()
        .map(|check| (check.text, check.settings))
        .collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_run_checks_partitions_need_an_alignment_file() {
    let config = config(
      AppCommand::Timetree,
      &json!({ "tree": "t.nwk", "metadata": "m.tsv", "sequence_length": 1000, "partition": [{ "name": "all" }] }),
    );
    assert_eq!(
      vec![
        "--partition needs --alignment: each entry takes sites of an --alignment file, but the run has no --alignment file"
      ],
      texts(&checks(AppCommand::Timetree, Some(&config), None, None))
    );
  }

  #[test]
  fn test_run_checks_partitions_check_the_model_of_a_run_without_alignments() {
    let config = config(
      AppCommand::Timetree,
      &json!({ "tree": "t.nwk", "metadata": "m.tsv", "sequence_length": 1000, "alphabet": "aa", "model": "jtt92" }),
    );
    assert_eq!(
      vec![(
        o!(
          "The model jtt92 needs the alphabet aa-no-stop because its matrix holds the 20 amino acids without the stop '*', but the alphabet is aa"
        ),
        vec_of_owned!["model", "model_params", "alphabet"],
      )],
      checks(AppCommand::Timetree, Some(&config), None, None)
        .into_iter()
        .map(|check| (check.text, check.settings))
        .collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_run_checks_partitions_wait_for_the_facts_of_every_alignment() {
    let config = config(
      AppCommand::Optimize,
      &json!({ "tree": "t.nwk", "alignment": ["a.fasta", "b.fasta"], "partition": [{ "alignment": 2 }] }),
    );
    let facts = InputFacts {
      alignments: vec![
        AlignmentFacts {
          min_length: 9,
          ..alignment(1, "a.fasta")
        },
        alignment(2, "b.fasta"),
      ],
      ..facts_with_tree(3)
    };
    let without_first = InputFacts {
      alignments: vec![alignment(2, "b.fasta")],
      ..facts_with_tree(3)
    };
    assert_eq!(
      (vec![o!("unaligned-sequences-1")], Vec::<String>::new()),
      (
        checks(AppCommand::Optimize, Some(&config), None, Some(&facts))
          .into_iter()
          .map(|check| check.id)
          .collect::<Vec<_>>(),
        checks(AppCommand::Optimize, Some(&config), None, Some(&without_first))
          .into_iter()
          .map(|check| check.id)
          .collect::<Vec<_>>()
      )
    );
  }

  #[test]
  fn test_run_checks_partitions_leave_an_entry_error_to_the_rejected_config() {
    let config = config(
      AppCommand::Ancestral,
      &json!({ "tree": "t.nwk", "alignment": ["a.fasta"], "partition": [{ "model_params": ["mu=2"] }] }),
    );
    let facts = InputFacts {
      alignments: vec![alignment(1, "a.fasta")],
      ..facts_with_tree(3)
    };
    let message = "partition 1: model_params needs a named model: infer fits the model from the data and takes no parameters, but got 'mu=2'";
    let settings = vec_of_owned!["partition"];
    let rejection = ConfigRejection {
      message,
      causes: &[],
      problems: &[],
      settings: &settings,
    };
    assert_eq!(
      vec![(o!("config-0"), o!(message))],
      checks(AppCommand::Ancestral, Some(&config), Some(rejection), Some(&facts))
        .into_iter()
        .map(|check| (check.id, check.text))
        .collect::<Vec<_>>()
    );
  }

  mod helpers {
    use crate::check_inputs::{AlignmentFacts, InputFacts, TreeFacts};
    use crate::check_inputs::{DateFacts, MetadataFacts};
    use crate::command::AppCommandExt;
    use crate::run_checks::{CheckContext, ConfigRejection, RunCheck, run_checks};
    use app_commands::command::AppCommand;
    use app_commands::config::source::ConfigProblem;
    use serde_json::{Map, Value};
    use std::path::PathBuf;
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
          duplicate_node_names: vec![],
        }),
        ..InputFacts::default()
      }
    }
    pub(super) fn alignment(partition: usize, path: &str) -> AlignmentFacts {
      AlignmentFacts {
        partition,
        path: PathBuf::from(path),
        sequences: 20,
        min_length: 10,
        max_length: 10,
        duplicate_names: vec![],
        tips_without_sequence: None,
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
