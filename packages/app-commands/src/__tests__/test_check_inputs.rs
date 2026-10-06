#[cfg(test)]
mod tests {
  use crate::check_inputs::{
    AlignmentFacts, DateFacts, InputKind, MetadataFacts, TreeFacts, check_inputs, metadata_summary,
  };
  use crate::command::AppCommand;
  use crate::config::properties::leaf_properties;
  use helpers::{data, request};
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
    let alignment = facts.alignment.clone().unwrap();
    assert_eq!(
      (
        Some(MetadataFacts {
          rows: 86,
          columns: vec!["name".to_owned(), "date".to_owned(), "country".to_owned()],
          id_column: "name".to_owned(),
          date_column: Some("date".to_owned()),
          dates: Some(DateFacts {
            readable: 86,
            unreadable: vec![],
            exact_days: 86,
            on_day_1_or_15: 54,
          }),
        }),
        AlignmentFacts {
          sequences: 86,
          min_length: alignment.max_length,
          max_length: alignment.max_length,
          duplicate_names: vec![],
        },
        Some(vec![]),
        Some(vec![]),
        0,
      ),
      (
        facts.metadata,
        alignment,
        facts.tips_without_metadata,
        facts.tips_without_sequence,
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
        }),
        (3, 4),
      ),
      (
        facts.tips_without_metadata,
        facts.tips_without_sequence,
        facts.tree.unwrap().polytomies,
        facts.metadata.unwrap().dates,
        (
          facts.alignment.as_ref().unwrap().min_length,
          facts.alignment.as_ref().unwrap().max_length
        ),
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
    assert_eq!(
      (true, None, None),
      (facts.tree.is_some(), facts.alignment, facts.tips_without_sequence)
    );
  }

  #[test]
  fn test_check_inputs_input_settings_are_settings_of_every_command_that_reads_them() {
    let missing = AppCommand::iter()
      .flat_map(|command| {
        let keys = leaf_properties(command.config_schema().as_value())
          .unwrap()
          .into_iter()
          .map(|leaf| leaf.key_path.join("."))
          .collect::<BTreeSet<_>>();
        let reads_metadata = command.inputs().iter().any(|input| input.kind == InputKind::Metadata);
        let metadata_keys: &[&str] = if reads_metadata {
          &["metadata_id_columns", "metadata_delimiters"]
        } else {
          &[]
        };
        let date_keys: &[&str] = if command.uses_dates() { &["date_column"] } else { &[] };
        command
          .inputs()
          .iter()
          .map(|input| input.kind.setting())
          .chain(metadata_keys.iter().copied())
          .chain(date_keys.iter().copied())
          .filter(|key| !keys.contains(*key))
          .map(|key| format!("{command}: {key}"))
          .collect_vec()
      })
      .collect_vec();
    assert_eq!(Vec::<String>::new(), missing);
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

    let command_dates = table().dates().unwrap();
    let read = metadata_summary(table());

    let date_facts = read.facts.dates.clone().unwrap();
    assert_eq!(
      (
        "strain".to_owned(),
        Some("date".to_owned()),
        command_dates.keys().cloned().collect::<BTreeSet<_>>(),
        command_dates.values().filter(|date| date.is_some()).count(),
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
    let read = metadata_summary(table);
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

  mod helpers {
    use crate::check_inputs::CheckInputsRequest;
    use crate::command::AppCommand;
    use crate::json_value::SparseConfig;
    use serde_json::Value;
    use std::path::{Path, PathBuf};

    pub(super) fn data(path: &str) -> PathBuf {
      Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data").join(path)
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
