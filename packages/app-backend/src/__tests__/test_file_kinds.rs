#[cfg(test)]
mod tests {
  use crate::file_kinds::{DetectedFile, FileKind, FileKindsRequest, file_kinds};
  use app_commands::command::AppCommand;
  use app_datasets::schema_directive;
  use helpers::{PackageWriter, kind_of, run_package, zip_with};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::fs;
  use tempfile::tempdir;
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::tree_named_as_alignment(     "tree.fasta",   "((A:1,B:1):1,C:1);\n",                       FileKind::Tree)]
  #[case::nexus_named_as_text(         "tree.txt",     "#NEXUS\nbegin trees;\ntree t = (A,B);\nend;\n", FileKind::Tree)]
  #[case::alignment_named_as_tree(     "aln.nwk",      ">A\nACGT\n>B\nACGA\n",                       FileKind::Alignment)]
  #[case::protein_alignment(           "proteins.txt", ">A\nMKVL\n>B\nMKEL\n",                       FileKind::Alignment)]
  #[case::metadata_named_as_tree(      "dates.nwk",    "strain\tdate\nA\t2020-01-01\n",              FileKind::Metadata)]
  #[case::comma_metadata_named_as_tsv( "dates.tsv",    "name,date\nA,2020.5\nB,2021.0\n",            FileKind::Metadata)]
  #[trace]
  fn test_file_kinds_reads_the_kind_from_the_content(#[case] name: &str, #[case] content: &str, #[case] expected: FileKind) {
    let dir = tempdir().unwrap();
    let path = dir.path().join(name);
    fs::write(&path, content).unwrap();
    assert_eq!(
      DetectedFile {
        path: path.clone(),
        kind: expected,
        commands: vec![],
        reason: None,
      },
      kind_of(&path)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::schema_directive(   "config.txt", Some("clock"), "tree: tree.nwk\n",                              vec![AppCommand::Clock])]
  #[case::keys_of_one_command("run.yaml",   None,          "method_anc: marginal\ntree: tree.nwk\n",         vec![AppCommand::Ancestral])]
  #[case::keys_of_two(        "a.json",     None,          "{\"tree\": \"tree.nwk\", \"clock_filter\": 4}", vec![AppCommand::Timetree, AppCommand::Clock])]
  #[trace]
  fn test_file_kinds_names_the_commands_of_a_configuration(
    #[case] name: &str,
    #[case] directive: Option<&str>,
    #[case] settings: &str,
    #[case] commands: Vec<AppCommand>,
  ) {
    let dir = tempdir().unwrap();
    let path = dir.path().join(name);
    let directive = directive.map_or_else(String::new, |command| format!("{}\n", schema_directive(command)));
    fs::write(&path, format!("{directive}{settings}")).unwrap();
    assert_eq!(
      DetectedFile {
        path: path.clone(),
        kind: FileKind::Configuration,
        commands,
        reason: None,
      },
      kind_of(&path)
    );
  }

  #[test]
  fn test_file_kinds_rejects_a_configuration_of_a_command_line_only_command() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("mugration.yaml");
    fs::write(
      &path,
      format!("{}\nattribute: country\n", schema_directive("mugration")),
    )
    .unwrap();
    assert_eq!(
      DetectedFile {
        path: path.clone(),
        kind: FileKind::Unknown,
        commands: vec![],
        reason: Some(o!(
          "The file name suggests a configuration, but no reader of TreeTime reads the file. Read as a configuration: the configuration is for `mugration`, which only the command-line tool runs; as an alignment: FASTA input is incorrectly formatted: expected at least one FASTA record starting with character '>', but none found; as a tree: line 1, column 3: expected end of input, comment, ')', ',', ':', ';', BEAST annotation; as a metadata table: Unable to detect metadata delimiter from candidates: ',', '\\t', ';'"
        )),
      },
      kind_of(&path)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::written_to_a_file(PackageWriter::File)]
  #[case::streamed(         PackageWriter::Stream)]
  #[trace]
  fn test_file_kinds_reads_a_run_package_of_the_package_writers_whatever_its_name(#[case] writer: PackageWriter) {
    let dir = tempdir().unwrap();
    let path = dir.path().join("download.bin");
    fs::write(&path, run_package(writer)).unwrap();
    assert_eq!(
      DetectedFile {
        path: path.clone(),
        kind: FileKind::RunPackage,
        commands: vec![],
        reason: None,
      },
      kind_of(&path)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zip_without_record(   &[("out/tree.nwk", "(A,B);")], "The file is not a run package: it holds no `run.json`")]
  #[case::zip_with_other_json(  &[("run.json", "{\"id\": 1}")],  "`run.json` of the run package is not a run record: When parsing JSON: Unexpected: unexpected unsigned integer, expected string at line 1 column 8 (path: id)")]
  #[case::record_in_two_folders(&[("a/run.json", "{}"), ("b/run.json", "{}")], "The file is not a run package: it holds no `run.json`")]
  #[trace]
  fn test_file_kinds_rejects_a_zip_archive_without_a_run_record(#[case] entries: &[(&str, &str)], #[case] reason: &str) {
    let dir = tempdir().unwrap();
    let path = dir.path().join("run.zip");
    fs::write(&path, zip_with(entries)).unwrap();
    assert_eq!(
      DetectedFile {
        path: path.clone(),
        kind: FileKind::Unknown,
        commands: vec![],
        reason: Some(o!(reason)),
      },
      kind_of(&path)
    );
  }

  #[test]
  fn test_file_kinds_gives_the_reason_of_each_reader_for_an_unknown_file() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("notes.nwk");
    fs::write(&path, "Samples arrive on Monday.\n").unwrap();
    assert_eq!(
      DetectedFile {
        path: path.clone(),
        kind: FileKind::Unknown,
        commands: vec![],
        reason: Some(o!(
          "The file name suggests a tree, but no reader of TreeTime reads the file. Read as a tree: line 1, column 9: expected end of input, comment, ')', ',', ':', ';', BEAST annotation; as a configuration: the document is not a mapping of settings; as an alignment: FASTA input is incorrectly formatted: expected at least one FASTA record starting with character '>', but none found; as a metadata table: Unable to detect metadata delimiter from candidates: ',', '\\t', ';'"
        )),
      },
      kind_of(&path)
    );
  }

  #[test]
  fn test_file_kinds_answers_each_file_in_the_order_of_the_request() {
    let dir = tempdir().unwrap();
    let tree = dir.path().join("a");
    let alignment = dir.path().join("b");
    let missing = dir.path().join("c");
    fs::write(&tree, "(A,B);").unwrap();
    fs::write(
      &alignment,
      indoc! {"
      >A
      ACGT
    "},
    )
    .unwrap();
    let kinds = file_kinds(&FileKindsRequest {
      files: vec![tree.clone(), alignment.clone(), missing.clone()],
    });
    assert_eq!(
      vec![
        (tree, FileKind::Tree),
        (alignment, FileKind::Alignment),
        (missing, FileKind::Unknown),
      ],
      kinds
        .files
        .into_iter()
        .map(|file| (file.path, file.kind))
        .collect::<Vec<_>>()
    );
  }

  mod helpers {
    use crate::__tests__::test_support::tests::{sparse, unconfined};
    use crate::file_kinds::{DetectedFile, FileKindsRequest, file_kinds};
    use crate::runs::manager::RunManager;
    use crate::runs::record::{CreateRunRequest, StartRunRequest};
    use app_commands::command::AppCommand;
    use serde_json::json;
    use std::io::{Cursor, Write};
    use std::path::Path;
    use tempfile::tempdir;
    use zip::ZipWriter;
    use zip::write::SimpleFileOptions;

    pub(super) fn kind_of(path: &Path) -> DetectedFile {
      let mut kinds = file_kinds(&FileKindsRequest {
        files: vec![path.to_path_buf()],
      });
      kinds.files.remove(0)
    }

    #[derive(Clone, Copy, Debug)]
    pub(super) enum PackageWriter {
      File,
      Stream,
    }

    pub(super) fn run_package(writer: PackageWriter) -> Vec<u8> {
      let root = tempdir().unwrap();
      let runs = RunManager::open(root.path()).unwrap();
      let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
      let created = runs
        .create(&CreateRunRequest {
          command: AppCommand::Clock,
          config: sparse(json!({ "tree": data.join("tree.nwk"), "metadata": data.join("metadata.tsv") })),
          defer_start: true,
          source_run_id: None,
          example: None,
        })
        .unwrap();
      runs
        .start(&created.id, &StartRunRequest::default(), unconfined())
        .unwrap()
        .run();
      let package = runs.package(&created.id).unwrap();
      match writer {
        PackageWriter::File => {
          let mut written = Cursor::new(vec![]);
          package.write(&mut written).unwrap();
          written.into_inner()
        },
        PackageWriter::Stream => {
          let mut streamed = vec![];
          package.stream(&mut streamed).unwrap();
          streamed
        },
      }
    }

    pub(super) fn zip_with(entries: &[(&str, &str)]) -> Vec<u8> {
      let mut writer = ZipWriter::new(Cursor::new(Vec::new()));
      for (name, content) in entries {
        writer.start_file(*name, SimpleFileOptions::default()).unwrap();
        writer.write_all(content.as_bytes()).unwrap();
      }
      writer.finish().unwrap().into_inner()
    }
  }
}
