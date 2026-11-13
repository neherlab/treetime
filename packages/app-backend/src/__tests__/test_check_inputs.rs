#[cfg(test)]
mod tests {
  use crate::check_inputs::{
    AlignmentFacts, DateFacts, DateSpan, InputKind, MetadataFacts, TreeFacts, check_inputs, metadata_summary,
  };
  use crate::command::AppCommandExt;
  use crate::input_detection::{
    DateForm, DateReading, DetectedAlphabet, DetectedDateFormat, DetectedDelimiter, DetectedDialect, DetectedLength,
    DetectedPartitions, DetectedSettings, PartitionFacts, SettingOrigin,
  };
  use crate::results::year_date::YearDate;
  use app_commands::command::AppCommand;
  use app_commands::commands::shared::alphabet::AlphabetNameCli;
  use app_commands::config::properties::setting_properties;
  use helpers::{data, detected, reading, request, span};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::json;
  use std::collections::BTreeSet;
  use std::fs;
  use strum::IntoEnumIterator;
  use tempfile::tempdir;
  use treetime_io::csv::{default_metadata_delimiters, default_name_candidates};
  use treetime_io::dates_csv::metadata_read;
  use treetime_utils::datetime::options::DateParserOptions;
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::iso(          DateForm::Calendar,  Some("%Y-%m-%d"),  "YYYY-MM-DD")]
  #[case::day_first(    DateForm::Calendar,  Some("%d/%m/%y"),  "DD/MM/YY")]
  #[case::month_name(   DateForm::Calendar,  Some("%b %e %Y"),  "MMM D YYYY")]
  #[case::literal(      DateForm::Calendar,  Some("%Y%%x"),     "YYYY%x")]
  #[case::no_format(    DateForm::Calendar,  None,              "calendar dates")]
  #[case::numeric(      DateForm::Numeric,   None,              "decimal years")]
  #[case::uncertain(    DateForm::Uncertain, None,              "uncertain dates")]
  #[case::range(        DateForm::Range,     None,              "date ranges")]
  #[trace]
  fn test_check_inputs_date_reading_text(#[case] form: DateForm, #[case] format: Option<&str>, #[case] expected: &str) {
    assert_eq!(expected, DateReading::new(form, format.map(str::to_owned), 1).text);
  }

  #[test]
  fn test_check_inputs_zika_86_tree_facts() {
    let facts = check_inputs(&request(
      AppCommand::Timetree,
      json!({ "tree": data("zika/86/tree.nwk") }),
    ))
    .unwrap();
    assert_eq!(
      Some(TreeFacts {
        tips: 86,
        internal_nodes: facts.tree.as_ref().unwrap().internal_nodes,
        polytomies: 15,
        unnamed_tips: 0,
        duplicate_node_names: vec![],
        tip_names: facts.tree.as_ref().unwrap().tip_names.clone(),
      }),
      facts.tree
    );
  }

  #[test]
  fn test_check_inputs_zika_86_metadata_and_alignment_facts() {
    let facts = check_inputs(&request(
      AppCommand::Timetree,
      json!({
        "tree": data("zika/86/tree.nwk"),
        "metadata": data("zika/86/metadata.tsv"),
        "alignment": [data("zika/86/aln.fasta.xz")],
      }),
    ))
    .unwrap();
    let alignment = facts.alignments[0].clone();
    let dates = facts.metadata.as_ref().unwrap().dates.clone().unwrap();
    assert_eq!(
      Some((o!("2013-10-25"), o!("2016-07-01"))),
      dates
        .span
        .as_ref()
        .map(|span| (span.earliest.date.clone().unwrap(), span.latest.date.clone().unwrap()))
    );
    assert_eq!(
      (
        Some(MetadataFacts {
          rows: 86,
          columns: vec![
            "name".to_owned(),
            "date".to_owned(),
            "country".to_owned(),
            "region".to_owned()
          ],
          id_column: "name".to_owned(),
          date_column: Some("date".to_owned()),
          dates: Some(DateFacts {
            readable: 86,
            unreadable: vec![],
            exact_days: 86,
            on_day_1_or_15: 54,
            span: dates.span,
          }),
        }),
        vec![AlignmentFacts {
          number: 1,
          path: data("zika/86/aln.fasta.xz"),
          sequences: 86,
          min_length: alignment.max_length,
          max_length: alignment.max_length,
          duplicate_names: vec![],
          tips_without_sequence: Some(vec![]),
        }],
        Some(vec![]),
        0,
      ),
      (
        facts.metadata,
        facts.alignments,
        facts.tips_without_metadata,
        facts.problems.len(),
      )
    );
  }

  #[test]
  fn test_check_inputs_reports_missing_samples_and_unreadable_dates() {
    let dir = tempdir().unwrap();
    fs::write(dir.path().join("tree.nwk"), "((A:1,B:1):1,(C:1,D:1,E:1):1);\n").unwrap();
    fs::write(
      dir.path().join("metadata.csv"),
      "strain,date\nA,2020-03-15\nB,2020-03-02\nC,soon\nX,2020-01-01\n",
    )
    .unwrap();
    fs::write(dir.path().join("aln.fasta"), ">A\nACGT\n>B\nACGT\n>C\nACG\n").unwrap();
    let facts = check_inputs(&request(
      AppCommand::Timetree,
      json!({
        "tree": dir.path().join("tree.nwk"),
        "metadata": dir.path().join("metadata.csv"),
        "alignment": [dir.path().join("aln.fasta")],
      }),
    ))
    .unwrap();
    assert_eq!(
      (
        Some(vec!["D".to_owned(), "E".to_owned()]),
        Some(vec!["D".to_owned(), "E".to_owned()]),
        1,
        Some(DateFacts {
          readable: 3,
          unreadable: vec!["C".to_owned()],
          exact_days: 3,
          on_day_1_or_15: 2,
          span: Some(DateSpan {
            earliest: YearDate::new(2020.0 + 0.5 / 366.0),
            latest: YearDate::new(2020.0 + 74.5 / 366.0),
          }),
        }),
        (3, 4),
      ),
      (
        facts.tips_without_metadata,
        facts.alignments[0].tips_without_sequence.clone(),
        facts.tree.unwrap().polytomies,
        facts.metadata.unwrap().dates,
        (facts.alignments[0].min_length, facts.alignments[0].max_length),
      )
    );
  }

  #[test]
  fn test_check_inputs_ignores_inputs_the_command_does_not_read() {
    let facts = check_inputs(&request(
      AppCommand::Mugration,
      json!({
        "tree": data("zika/86/tree.nwk"),
        "alignment": [data("zika/86/aln.fasta.xz")],
      }),
    ))
    .unwrap();
    assert_eq!((true, vec![]), (facts.tree.is_some(), facts.alignments));
  }

  #[test]
  fn test_check_inputs_reports_facts_per_alignment_file() {
    let dir = tempdir().unwrap();
    fs::create_dir_all(dir.path().join("ha")).unwrap();
    fs::create_dir_all(dir.path().join("na")).unwrap();
    let tree = dir.path().join("tree.nwk");
    let ha = dir.path().join("ha/aligned.fasta");
    let na = dir.path().join("na/aligned.fasta");
    let missing = dir.path().join("mp/aligned.fasta");
    fs::write(&tree, "((A:1,B:1):1,C:1);\n").unwrap();
    fs::write(&ha, ">A\nACGT\n>A\nACGA\n>B\nACGT\n").unwrap();
    fs::write(&na, ">A\nAC\n>C\nACG\n").unwrap();

    let facts = check_inputs(&request(
      AppCommand::Optimize,
      json!({ "tree": tree, "alignment": [ha, na, missing] }),
    ))
    .unwrap();

    assert_eq!(
      (
        vec![
          AlignmentFacts {
            number: 1,
            path: ha,
            sequences: 3,
            min_length: 4,
            max_length: 4,
            duplicate_names: vec!["A".to_owned()],
            tips_without_sequence: Some(vec!["C".to_owned()]),
          },
          AlignmentFacts {
            number: 2,
            path: na,
            sequences: 2,
            min_length: 2,
            max_length: 3,
            duplicate_names: vec![],
            tips_without_sequence: Some(vec!["B".to_owned()]),
          },
        ],
        vec![(InputKind::Alignment, Some(3))],
      ),
      (
        facts.alignments,
        facts
          .problems
          .iter()
          .map(|problem| (problem.input, problem.alignment))
          .collect::<Vec<_>>()
      )
    );
  }

  #[test]
  fn test_check_inputs_sequence_lengths_count_the_first_record_of_each_name() {
    let dir = tempdir().unwrap();
    let tree = dir.path().join("tree.nwk");
    let aln = dir.path().join("aln.fasta");
    fs::write(&tree, "(A:1,B:1);\n").unwrap();
    fs::write(&aln, ">A\nACGT\n>A\nAC\n>B\nACGT\n").unwrap();

    let facts = check_inputs(&request(
      AppCommand::Optimize,
      json!({ "tree": tree, "alignment": [aln] }),
    ))
    .unwrap();

    assert_eq!(
      vec![(4, 4, vec!["A".to_owned()])],
      facts
        .alignments
        .iter()
        .map(|file| (file.min_length, file.max_length, file.duplicate_names.clone()))
        .collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_check_inputs_input_settings_are_settings_of_every_command_that_reads_them() {
    let missing = AppCommand::iter()
      .flat_map(|command| {
        let keys = setting_properties(command.config_schema().as_value())
          .unwrap()
          .into_iter()
          .map(|setting| setting.key)
          .collect::<BTreeSet<_>>();
        let reads_metadata = command.inputs().iter().any(|input| input.kind == InputKind::Metadata);
        let metadata_keys: &[&str] = if reads_metadata {
          &["metadata_id_columns", "metadata_delimiters"]
        } else {
          &[]
        };
        let date_keys: &[&str] = if command.uses_dates() {
          &["date_column", "date_format"]
        } else {
          &[]
        };
        let reads_tree = command.inputs().iter().any(|input| input.kind == InputKind::Tree);
        let tree_keys: &[&str] = if reads_tree { &["newick_dialect"] } else { &[] };
        command
          .inputs()
          .iter()
          .map(|input| input.kind.setting())
          .chain(metadata_keys.iter().copied())
          .chain(date_keys.iter().copied())
          .chain(tree_keys.iter().copied())
          .filter(|key| !keys.contains(*key))
          .map(|key| format!("{command}: {key}"))
          .collect_vec()
      })
      .collect_vec();
    assert_eq!(Vec::<String>::new(), missing);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::network_dialect(  json!("enewick,beast"), (None, vec![InputKind::Tree]))]
  #[case::default_dialect(  json!(null),            (Some(4), vec![]))]
  #[trace]
  fn test_check_inputs_reads_the_tree_in_the_configured_dialect(
    #[case] newick_dialect: serde_json::Value,
    #[case] expected: (Option<usize>, Vec<InputKind>),
  ) {
    let dir = tempdir().unwrap();
    let tree = dir.path().join("tree.nwk");
    fs::write(&tree, "((A:1,(B:1)#H0:0.5):1,(#H0:0.7,C:1):1);\n").unwrap();
    let mut config = json!({ "tree": tree });
    if !newick_dialect.is_null() {
      config["newick_dialect"] = newick_dialect;
    }

    let facts = check_inputs(&request(AppCommand::Clock, config)).unwrap();

    assert_eq!(
      expected,
      (
        facts.tree.map(|tree| tree.tips),
        facts.problems.iter().map(|problem| problem.input).collect::<Vec<_>>()
      )
    );
  }

  #[test]
  fn test_check_inputs_reports_unreadable_inputs_as_problems() {
    let dir = tempdir().unwrap();
    let facts = check_inputs(&request(
      AppCommand::Timetree,
      json!({ "tree": dir.path().join("missing.nwk") }),
    ))
    .unwrap();
    assert_eq!(
      (None, vec![InputKind::Tree]),
      (
        facts.tree,
        facts.problems.iter().map(|problem| problem.input).collect::<Vec<_>>()
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::csv(             b',',  ',')]
  #[case::tsv(             b'\t', '\t')]
  #[case::tsv_named_as_csv(b',',  '\t')]
  #[case::csv_named_as_tsv(b'\t', ',')]
  #[trace]
  fn test_check_inputs_metadata_facts_match_the_command_dates(#[case] path_delimiter: u8, #[case] delimiter: char) {
    let rows = [["strain", "date", "country"], ["A", "2020-01-15", "usa"], ["B", "2020-02-03", "peru"], ["C", "soon", "chile"]];
    let text = rows.iter().map(|row| row.join(&delimiter.to_string())).join("\n");
    let table = || {
      metadata_read(
        text.as_bytes(),
        Some(path_delimiter),
        &default_metadata_delimiters(),
        &default_name_candidates(),
        None,
        None,
      )
      .unwrap()
    };

    let command_dates = table().dates(&DateParserOptions::default()).unwrap();
    let read = metadata_summary(table(), &DateParserOptions::default());

    let date_facts = read.facts.dates.clone().unwrap();
    assert_eq!(
      (
        "strain".to_owned(),
        Some("date".to_owned()),
        command_dates.iter().map(|(name, _)| name.clone()).collect::<BTreeSet<_>>(),
        command_dates.iter().filter(|(_, date)| date.is_some()).count(),
        vec!["C".to_owned()],
        1,
      ),
      (
        read.facts.id_column,
        read.facts.date_column,
        read.names,
        date_facts.readable,
        date_facts.unreadable,
        date_facts.on_day_1_or_15,
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::without_format(json!(null),         DateFacts { readable: 0, unreadable: vec!["A".to_owned()], exact_days: 0, on_day_1_or_15: 0, span: None })]
  #[case::with_format(   json!("%d %b %Y"),   DateFacts { readable: 1, unreadable: vec![],               exact_days: 1, on_day_1_or_15: 1, span: Some(span(2016.0 + 166.5 / 366.0)) })]
  #[trace]
  fn test_check_inputs_reads_dates_in_the_configured_date_format(
    #[case] date_format: serde_json::Value,
    #[case] expected: DateFacts,
  ) {
    let dir = tempdir().unwrap();
    let path = dir.path().join("metadata.tsv");
    fs::write(&path, "strain\tdate\nA\t15 Jun 2016\n").unwrap();

    let facts = check_inputs(&request(
      AppCommand::Clock,
      json!({ "metadata": path, "date_format": date_format }),
    ))
    .unwrap();

    assert_eq!(Some(expected), facts.metadata.and_then(|metadata| metadata.dates));
  }

  #[test]
  fn test_check_inputs_metadata_without_date_column_has_no_date_facts() {
    let table = metadata_read(
      &b"name\tcountry\nA\tusa\n"[..],
      Some(b'\t'),
      &default_metadata_delimiters(),
      &default_name_candidates(),
      None,
      None,
    )
    .unwrap();
    let read = metadata_summary(table, &DateParserOptions::default());
    assert_eq!(
      (1, "name".to_owned(), None, None, true),
      (
        read.facts.rows,
        read.facts.id_column,
        read.facts.date_column,
        read.facts.dates,
        read.date_problem.is_none()
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::iso_dates(              "zika/20",     ("name",      vec![reading(DateForm::Calendar,  Some("%Y-%m-%d"), 20)],                                                  20, 20), 10807)]
  #[case::numeric_dates(          "ebola/20",    ("name",      vec![reading(DateForm::Numeric,   None,             20)],                                                  20, 20), 19006)]
  #[case::numeric_dates_missing(  "flu/h3n2/20", ("name",      vec![reading(DateForm::Numeric,   None,             19)],                                                  19, 19), 1407)]
  #[case::iso_and_uncertain(      "lassa/L/20",  ("accession", vec![reading(DateForm::Calendar,  Some("%Y-%m-%d"), 19), reading(DateForm::Uncertain, None, 1)],           20, 20), 7889)]
  #[case::uncertain_dates(        "tb/20",       ("strain",    vec![reading(DateForm::Uncertain, None,             20)],                                                  20, 20), 216)]
  #[case::mostly_uncertain(       "dengue/20",   ("strain",    vec![reading(DateForm::Uncertain, None,             16), reading(DateForm::Calendar, Some("%Y-%m-%d"), 4)],  20, 20), 10723)]
  #[trace]
  fn test_check_inputs_detects_the_data_settings_of_example_datasets(
    #[case] dataset: &str,
    #[case] (id_column, readings, read, rows): (&str, Vec<DateReading>, usize, usize),
    #[case] length: usize,
  ) {
    let facts = check_inputs(&request(
      AppCommand::Timetree,
      json!({
        "tree": data(&format!("{dataset}/tree.nwk")),
        "metadata": data(&format!("{dataset}/metadata.tsv")),
        "alignment": [data(&format!("{dataset}/aln.fasta.xz"))],
      }),
    ))
    .unwrap();
    assert_eq!(
      DetectedSettings {
        newick_dialect: Some(DetectedDialect {
          dialect: o!("classic,beast"),
          origin: SettingOrigin::Default,
          setting: o!("newick_dialect"),
        }),
        metadata_delimiter: Some(DetectedDelimiter {
          delimiter: o!("\t"),
          origin: SettingOrigin::Detected,
          setting: o!("metadata_delimiters"),
        }),
        id_column: Some(detected(id_column, SettingOrigin::Detected, "metadata_id_columns")),
        date_column: Some(detected("date", SettingOrigin::Detected, "date_column")),
        date_format: Some(DetectedDateFormat {
          readings,
          read,
          rows,
          origin: SettingOrigin::Detected,
          setting: o!("date_format"),
        }),
        alphabet: Some(DetectedAlphabet {
          alphabet: AlphabetNameCli::Nuc,
          origin: SettingOrigin::Default,
          setting: o!("alphabet"),
        }),
        sequence_length: Some(DetectedLength {
          length,
          origin: SettingOrigin::Detected,
          setting: o!("alignment"),
        }),
        partitions: Some(DetectedPartitions {
          partitions: vec![PartitionFacts {
            label: o!("partition 1"),
            alignment: 1,
            sites: length,
            model: o!("infer"),
          }],
          origin: SettingOrigin::Default,
          setting: o!("partition"),
        }),
      },
      facts.detected
    );
  }

  #[test]
  fn test_check_inputs_detected_settings_name_the_settings_that_set_them() {
    let dir = tempdir().unwrap();
    let tree = dir.path().join("tree.nwk");
    let metadata = dir.path().join("metadata.txt");
    let aln = dir.path().join("aln.fasta");
    fs::write(&tree, "((A:1,B:1):1,C:1);\n").unwrap();
    fs::write(
      &metadata,
      "id;country;sampled\nA;usa;15 Jun 2016\nB;peru;2016.5\nC;chile;2016-06-15\n",
    )
    .unwrap();
    fs::write(&aln, ">A\nMKVL\n>B\nMKVI\n>C\nMRVL\n").unwrap();

    let facts = check_inputs(&request(
      AppCommand::Timetree,
      json!({
        "tree": tree,
        "metadata": metadata,
        "alignment": [aln],
        "newick_dialect": "classic,nhx",
        "metadata_delimiters": [";"],
        "metadata_id_columns": ["id"],
        "date_column": "sampled",
        "date_format": "%d %b %Y",
        "alphabet": "aa",
        "partition": [{ "name": "start", "sites": ["1-2"] }, { "name": "end" }],
      }),
    ))
    .unwrap();

    assert_eq!(
      (
        vec![],
        DetectedSettings {
          newick_dialect: Some(DetectedDialect {
            dialect: o!("classic,nhx"),
            origin: SettingOrigin::Setting,
            setting: o!("newick_dialect"),
          }),
          metadata_delimiter: Some(DetectedDelimiter {
            delimiter: o!(";"),
            origin: SettingOrigin::Setting,
            setting: o!("metadata_delimiters"),
          }),
          id_column: Some(detected("id", SettingOrigin::Setting, "metadata_id_columns")),
          date_column: Some(detected("sampled", SettingOrigin::Setting, "date_column")),
          date_format: Some(DetectedDateFormat {
            readings: vec![
              reading(DateForm::Numeric, None, 1),
              reading(DateForm::Calendar, Some("%Y-%m-%d"), 1),
              reading(DateForm::Calendar, Some("%d %b %Y"), 1),
            ],
            read: 3,
            rows: 3,
            origin: SettingOrigin::Setting,
            setting: o!("date_format"),
          }),
          alphabet: Some(DetectedAlphabet {
            alphabet: AlphabetNameCli::Aa,
            origin: SettingOrigin::Setting,
            setting: o!("alphabet"),
          }),
          sequence_length: Some(DetectedLength {
            length: 4,
            origin: SettingOrigin::Detected,
            setting: o!("alignment"),
          }),
          partitions: Some(DetectedPartitions {
            partitions: vec![
              PartitionFacts {
                label: o!("partition 1 'start'"),
                alignment: 1,
                sites: 2,
                model: o!("infer"),
              },
              PartitionFacts {
                label: o!("partition 2 'end'"),
                alignment: 1,
                sites: 2,
                model: o!("infer"),
              },
            ],
            origin: SettingOrigin::Setting,
            setting: o!("partition"),
          }),
        },
      ),
      (
        facts
          .problems
          .iter()
          .map(|problem| problem.message.clone())
          .collect::<Vec<_>>(),
        facts.detected
      )
    );
  }

  #[test]
  fn test_check_inputs_sequence_length_without_alignment_comes_from_its_setting() {
    let facts = check_inputs(&request(
      AppCommand::Timetree,
      json!({ "tree": data("zika/20/tree.nwk"), "sequence_length": 10807 }),
    ))
    .unwrap();
    assert_eq!(
      (
        Some(DetectedLength {
          length: 10807,
          origin: SettingOrigin::Setting,
          setting: o!("sequence_length"),
        }),
        None,
        None,
      ),
      (
        facts.detected.sequence_length,
        facts.detected.alphabet,
        facts.detected.partitions
      )
    );
  }

  #[test]
  fn test_check_inputs_detects_only_the_settings_of_the_command() {
    let facts = check_inputs(&request(
      AppCommand::Clock,
      json!({ "tree": data("zika/20/tree.nwk"), "metadata": data("zika/20/metadata.tsv") }),
    ))
    .unwrap();
    assert_eq!(
      (true, true, true, None, None, None),
      (
        facts.detected.newick_dialect.is_some(),
        facts.detected.id_column.is_some(),
        facts.detected.date_format.is_some(),
        facts.detected.alphabet,
        facts.detected.sequence_length,
        facts.detected.partitions,
      )
    );
  }

  #[test]
  fn test_check_inputs_unaligned_sequences_have_no_sequence_length() {
    let dir = tempdir().unwrap();
    let aln = dir.path().join("aln.fasta");
    fs::write(&aln, ">A\nACGT\n>B\nACG\n").unwrap();
    let facts = check_inputs(&request(AppCommand::Ancestral, json!({ "alignment": [aln] }))).unwrap();
    assert_eq!(
      (None, None),
      (facts.detected.sequence_length, facts.detected.partitions)
    );
  }

  mod helpers {
    use crate::check_inputs::{CheckInputsRequest, DateSpan};
    use crate::input_detection::{DateForm, DateReading, DetectedColumn, SettingOrigin};
    use crate::results::year_date::YearDate;
    use app_commands::command::AppCommand;
    use app_commands::json_value::SparseConfig;
    use serde_json::Value;
    use std::path::{Path, PathBuf};

    pub(super) fn data(path: &str) -> PathBuf {
      Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data").join(path)
    }

    pub(super) fn span(year: f64) -> DateSpan {
      DateSpan {
        earliest: YearDate::new(year),
        latest: YearDate::new(year),
      }
    }

    pub(super) fn reading(form: DateForm, format: Option<&str>, count: usize) -> DateReading {
      DateReading::new(form, format.map(str::to_owned), count)
    }

    pub(super) fn detected(column: &str, origin: SettingOrigin, setting: &str) -> DetectedColumn {
      DetectedColumn {
        column: column.to_owned(),
        origin,
        setting: setting.to_owned(),
      }
    }

    pub(super) fn request(command: AppCommand, config: Value) -> CheckInputsRequest {
      let Value::Object(config) = config else {
        panic!("a test configuration must be a mapping");
      };
      CheckInputsRequest {
        command,
        config: SparseConfig(config),
      }
    }
  }
}
