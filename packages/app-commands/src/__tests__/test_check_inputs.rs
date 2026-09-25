#[cfg(test)]
mod tests {
  use crate::check_inputs::{
    AlignmentFacts, CheckInputsRequest, DateFacts, InputKind, MetadataFacts, TreeFacts, check_inputs, metadata_read,
  };
  use helpers::data;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use std::fs;
  use tempfile::tempdir;
  use treetime_io::csv::{default_metadata_delimiters, default_name_candidates};
  use treetime_io::dates_csv::read_metadata_table_from_reader;

  #[test]
  fn test_check_inputs_zika_86_tree_facts() {
    let facts = check_inputs(&CheckInputsRequest {
      tree: Some(data("zika/86/tree.nwk")),
      ..CheckInputsRequest::default()
    });
    assert_eq!(
      Some(TreeFacts {
        tips: 86,
        internal_nodes: facts.tree.as_ref().unwrap().internal_nodes,
        polytomies: 15,
        unnamed_tips: 0,
        duplicate_tip_names: vec![],
      }),
      facts.tree
    );
  }

  #[test]
  fn test_check_inputs_zika_86_metadata_and_alignment_facts() {
    let facts = check_inputs(&CheckInputsRequest {
      tree: Some(data("zika/86/tree.nwk")),
      metadata: Some(data("zika/86/metadata.tsv")),
      alignment: vec![data("zika/86/aln.fasta.xz")],
      ..CheckInputsRequest::default()
    });
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
    let facts = check_inputs(&CheckInputsRequest {
      tree: Some(dir.path().join("tree.nwk")),
      metadata: Some(dir.path().join("metadata.csv")),
      alignment: vec![dir.path().join("aln.fasta")],
      ..CheckInputsRequest::default()
    });
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
  fn test_check_inputs_reports_unreadable_inputs_as_problems() {
    let dir = tempdir().unwrap();
    let facts = check_inputs(&CheckInputsRequest {
      tree: Some(dir.path().join("missing.nwk")),
      ..CheckInputsRequest::default()
    });
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
  #[case::csv(             "metadata.csv", ',')]
  #[case::tsv(             "metadata.tsv", '\t')]
  #[case::tsv_named_as_csv("metadata.csv", '\t')]
  #[case::csv_named_as_tsv("metadata.tsv", ',')]
  #[trace]
  fn test_check_inputs_metadata_facts_match_the_command_dates(#[case] file_name: &str, #[case] delimiter: char) {
    let rows = [["strain", "date", "country"], ["A", "2020-01-15", "usa"], ["B", "2020-02-03", "peru"], ["C", "soon", "chile"]];
    let text = rows.iter().map(|row| row.join(&delimiter.to_string())).join("\n");
    let table = || {
      read_metadata_table_from_reader(
        text.as_bytes(),
        file_name,
        &default_metadata_delimiters(),
        &default_name_candidates(),
        None,
        None,
      )
      .unwrap()
    };

    let command_dates = table().dates().unwrap();
    let read = metadata_read(table());

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
    let table = read_metadata_table_from_reader(
      &b"name\tcountry\nA\tusa\n"[..],
      "traits.tsv",
      &default_metadata_delimiters(),
      &default_name_candidates(),
      None,
      None,
    )
    .unwrap();
    let read = metadata_read(table);
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
    use std::path::{Path, PathBuf};

    pub(super) fn data(path: &str) -> PathBuf {
      Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data").join(path)
    }
  }
}
