#[cfg(test)]
mod tests {
  use crate::check_inputs::{
    AlignmentFacts, CheckInputsRequest, DateFacts, InputKind, MetadataFacts, TreeFacts, check_inputs,
  };
  use helpers::data;
  use pretty_assertions::assert_eq;
  use std::fs;
  use tempfile::tempdir;

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
          id_column: Some("name".to_owned()),
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

  #[test]
  fn test_check_inputs_metadata_without_date_column_has_no_date_facts() {
    let dir = tempdir().unwrap();
    fs::write(dir.path().join("traits.tsv"), "name\tcountry\nA\tusa\n").unwrap();
    let facts = check_inputs(&CheckInputsRequest {
      metadata: Some(dir.path().join("traits.tsv")),
      ..CheckInputsRequest::default()
    });
    let metadata = facts.metadata.unwrap();
    assert_eq!(
      (1, Some("name".to_owned()), None, None, 0),
      (
        metadata.rows,
        metadata.id_column,
        metadata.date_column,
        metadata.dates,
        facts.problems.len()
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
