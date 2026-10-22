#[cfg(test)]
mod tests {
  use helpers::{ALIGNMENT, Run, TREE, path_in, planned_lines, treetime_in, write_file, yaml_block};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use tempfile::tempdir;

  #[test]
  fn test_cli_dry_run_prints_a_config_that_runs_from_another_folder() {
    let inputs = tempdir().unwrap();
    write_file(inputs.path(), "tree.nwk", TREE);
    write_file(inputs.path(), "aln.fasta", ALIGNMENT);
    let dry_run = treetime_in(
      inputs.path(),
      &[
        "ancestral",
        "--tree",
        "tree.nwk",
        "--alignment",
        "aln.fasta",
        "--output-all",
        "out",
        "--dry-run",
      ],
    );
    assert_eq!((Some(0), String::new()), (dry_run.code, dry_run.stderr.clone()));
    assert!(!inputs.path().join("out").exists());

    let elsewhere = tempdir().unwrap();
    write_file(elsewhere.path(), "run.yaml", &dry_run.stdout);
    let again = treetime_in(elsewhere.path(), &["ancestral", "--config", "run.yaml", "--dry-run"]);
    assert_eq!(dry_run, again);

    let run = treetime_in(elsewhere.path(), &["ancestral", "--config", "run.yaml"]);
    assert_eq!(
      (Some(0), true),
      (run.code, inputs.path().join("out").join("ancestral.nwk").is_file())
    );
  }

  #[test]
  fn test_cli_dry_run_with_config_and_jobs_reads_no_input() {
    let dir = tempdir().unwrap();
    write_file(
      dir.path(),
      "run.yaml",
      indoc! {"
        tree: missing.nwk
        alignment:
          - missing.fasta
        output_all: out
      "},
    );
    let run = treetime_in(
      dir.path(),
      &["ancestral", "--config", "run.yaml", "--dry-run", "-j", "1"],
    );
    let tree_line = format!("tree: \"{}\"", path_in(dir.path(), "missing.nwk"));
    let nwk_line = format!("# {}", path_in(dir.path(), "out/ancestral.nwk"));
    assert_eq!(
      (Some(0), String::new(), true, true, false),
      (
        run.code,
        run.stderr,
        run.stdout.lines().any(|line| line == tree_line),
        run.stdout.lines().any(|line| line == nwk_line),
        dir.path().join("out").exists(),
      )
    );
  }

  #[test]
  fn test_cli_dry_run_lists_the_conditional_output_of_auto_branch_lengths() {
    let dir = tempdir().unwrap();
    let run = treetime_in(
      dir.path(),
      &[
        "timetree",
        "--tree",
        "tree.nwk",
        "--metadata",
        "metadata.tsv",
        "--alignment",
        "aln.fasta",
        "--output-all",
        "out",
        "--dry-run",
      ],
    );
    let comments = run
      .stdout
      .lines()
      .filter(|line| line.starts_with('#'))
      .filter(|line| !line.starts_with("# /") || line.contains("sequences"))
      .map(str::to_owned)
      .collect::<Vec<_>>();
    assert_eq!(
      (
        Some(0),
        true,
        vec![
          "# branch_length_mode: auto resolves to marginal or input from the input tree when the run starts; the \
           settings that need re-estimated branch lengths are checked then"
            .to_owned(),
          "# Planned output files:".to_owned(),
          format!(
            "# {} (written when auto resolves to marginal)",
            path_in(dir.path(), "out/timetree.sequences.fasta")
          ),
        ]
      ),
      (
        run.code,
        run.stdout.lines().any(|line| line == "branch_length_mode: \"auto\""),
        comments
      )
    );
  }

  #[test]
  fn test_cli_dry_run_prints_partition_entries_as_written() {
    let dir = tempdir().unwrap();
    let dry_run = treetime_in(
      dir.path(),
      &[
        "ancestral",
        "--tree",
        "missing.nwk",
        "--alignment",
        "missing.fasta",
        "--partition",
        r"{name: cp3, sites: ['3-1500\3'], model: hky85, model_params: ['kappa=4']}",
        "--partition",
        "{name: rest}",
        "--output-all",
        "out",
        "--dry-run",
      ],
    );
    let expected = indoc! {r#"
      partition:
        - name: "cp3"
          sites:
            - "3-1500\\3"
          model: "hky85"
          model_params:
            - "kappa=4"
        - name: "rest"
    "#};
    write_file(dir.path(), "run.yaml", &dry_run.stdout);
    let again = treetime_in(dir.path(), &["ancestral", "--config", "run.yaml", "--dry-run"]);
    assert_eq!(
      (Some(0), String::new(), expected.to_owned(), dry_run.stdout.clone()),
      (
        dry_run.code,
        dry_run.stderr,
        yaml_block(&dry_run.stdout, "partition"),
        again.stdout
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::names_from_annotation(&["--annotation", "genes.gff"],          &["aa/{cds}.fasta"])]
  #[case::names_from_selection( &["--cds-selection", "S", "N"],          &["aa/N.fasta", "aa/S.fasta"])]
  #[trace]
  fn test_cli_dry_run_lists_translations_by_where_the_cds_names_come_from(
    #[case] cds_args: &[&str],
    #[case] expected_translations: &[&str],
  ) {
    let dir = tempdir().unwrap();
    let args = [
      &["ancestral", "--tree", "tree.nwk", "--alignment", "aln.fasta", "--translations", "tr/{cds}.fasta"],
      cds_args,
      &["--output-sequences", "-", "--output-translations", "aa/{cds}.fasta", "--dry-run"],
    ]
    .concat();
    let run = treetime_in(dir.path(), &args);
    let expected = ["# Planned output files:".to_owned(), "# -".to_owned()]
      .into_iter()
      .chain(expected_translations.iter().map(|path| format!("# {}", path_in(dir.path(), path))))
      .collect::<Vec<_>>();
    assert_eq!((Some(0), expected), (run.code, planned_lines(&run.stdout)));
  }

  #[test]
  fn test_cli_dry_run_writes_each_planned_path_on_one_line() {
    let dir = tempdir().unwrap();
    let dry_run = treetime_in(
      dir.path(),
      &[
        "ancestral",
        "--tree",
        "tree.nwk",
        "--alignment",
        "aln.fasta",
        "--output-tree-nwk",
        "out\nnext.nwk",
        "--dry-run",
      ],
    );
    write_file(dir.path(), "run.yaml", &dry_run.stdout);
    let again = treetime_in(dir.path(), &["ancestral", "--config", "run.yaml", "--dry-run"]);
    assert_eq!(
      (
        Some(0),
        vec![
          "# Planned output files:".to_owned(),
          format!("# {}", path_in(dir.path(), "out\\nnext.nwk"))
        ],
        dry_run.stdout.clone()
      ),
      (again.code, planned_lines(&dry_run.stdout), again.stdout)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::prune_collapse_empty_without_alignment(
    &["prune", "--tree", "tree.nwk", "--collapse-empty", "--output-all", "out"],
    "The --collapse-empty requires --alignment. Without sequence data, it's not possible to determine which branches \
     lack mutations.")]
  #[case::ancestral_without_alignment(
    &["ancestral", "--tree", "tree.nwk", "--output-all", "out"],
    "--alignment is required: pass one or more FASTA files, or '-' to read standard input")]
  #[case::ancestral_standard_input_twice(
    &["ancestral", "--tree", "tree.nwk", "--alignment", "-", "--alignment", "-", "--output-all", "out"],
    "Standard input can be read only once, but --alignment lists it 2 times")]
  #[case::ancestral_annotation_without_translations(
    &["ancestral", "--tree", "tree.nwk", "--alignment", "aln.fasta", "--annotation", "genes.gff", "--output-all", "out"],
    "--translations is required when using --cds-selection, --annotation, or --aa-root-sequence")]
  #[case::ancestral_translations_with_two_alignments(
    &["ancestral", "--tree", "tree.nwk", "--alignment", "a.fasta", "--alignment", "b.fasta", "--translations", "tr/{cds}.fasta",
      "--cds-selection", "S", "--output-all", "out"],
    "--translations and --annotation need exactly one --alignment file, but 2 were given: CDS coordinates refer to \
     the positions of one nucleotide alignment")]
  #[case::ancestral_translation_output_without_translations(
    &["ancestral", "--tree", "tree.nwk", "--alignment", "aln.fasta", "--output-translations", "aa/{cds}.fasta"],
    "--output-translations requires --translations, which names the amino-acid alignments to reconstruct")]
  #[case::ancestral_translation_output_without_placeholder(
    &["ancestral", "--tree", "tree.nwk", "--alignment", "aln.fasta", "--translations", "tr/{cds}.fasta",
      "--cds-selection", "S", "N", "--output-translations", "aa.fasta"],
    "--output-translations needs the placeholder '{cds}' or '%GENE' in its file name when reconstructing several \
     CDSes, otherwise the output of each CDS overwrites the same file")]
  #[case::timetree_standard_input_twice(
    &["timetree", "--tree", "tree.nwk", "--metadata", "dates.tsv", "--alignment", "-", "--alignment", "-",
      "--output-all", "out"],
    "Standard input can be read only once, but --alignment lists it 2 times")]
  #[case::prune_standard_input_twice(
    &["prune", "--tree", "tree.nwk", "--alignment", "-", "--alignment", "-", "--output-all", "out"],
    "Standard input can be read only once, but --alignment lists it 2 times")]
  #[case::optimize_without_alignment(
    &["optimize", "--tree", "tree.nwk", "--output-all", "out"],
    "--alignment is required: pass one or more FASTA files, or '-' to read standard input")]
  #[trace]
  fn test_cli_dry_run_reports_the_argument_errors_of_a_run(#[case] args: &[&str], #[case] message: &str) {
    let dir = tempdir().unwrap();
    let expected = Run { code: Some(2), stdout: String::new(), stderr: format!("error: {message}\n") };
    let run = treetime_in(dir.path(), args);
    let dry_run = treetime_in(dir.path(), &[args, &["--dry-run"]].concat());
    assert_eq!((&expected, &expected), (&run, &dry_run));
  }

  #[test]
  fn test_cli_dry_run_stops_on_an_argument_error() {
    let dir = tempdir().unwrap();
    let run = treetime_in(
      dir.path(),
      &[
        "timetree",
        "--tree",
        "tree.nwk",
        "--metadata",
        "metadata.tsv",
        "--branch-length-mode",
        "marginal",
        "--output-all",
        "out",
        "--dry-run",
      ],
    );
    assert_eq!(
      Run {
        code: Some(2),
        stdout: String::new(),
        stderr: "error: --branch-length-mode=marginal requires an alignment: add --alignment or use \
                 --branch-length-mode=input\n"
          .to_owned(),
      },
      run
    );
  }

  #[test]
  fn test_cli_examples_get_dry_run_prints_the_download_without_downloading() {
    let dir = tempdir().unwrap();
    let run = treetime_in(
      dir.path(),
      &[
        "examples",
        "get",
        "--output-dir",
        "examples",
        "--url",
        "https://example.org/examples.zip",
        "--dry-run",
      ],
    );
    assert_eq!(
      (
        Run {
          code: Some(0),
          stdout: format!(
            "Download: https://example.org/examples.zip\nUnpack into: {}\n",
            path_in(dir.path(), "examples")
          ),
          stderr: String::new(),
        },
        false
      ),
      (run, dir.path().join("examples").exists())
    );
  }

  #[test]
  fn test_cli_examples_get_dry_run_rejects_a_folder_that_is_not_empty() {
    let dir = tempdir().unwrap();
    write_file(dir.path(), "existing.txt", "");
    let folder = dir.path().to_str().unwrap();
    let run = treetime_in(dir.path(), &["examples", "get", "--output-dir", folder, "--dry-run"]);
    let message = format!("the folder '{folder}' is not empty");
    assert_eq!(
      (Some(1), String::new(), true),
      (run.code, run.stdout, run.stderr.contains(&message))
    );
  }

  #[test]
  fn test_cli_schema_write_dry_run_lists_the_files_without_writing() {
    let dir = tempdir().unwrap();
    let run = treetime_in(
      dir.path(),
      &["schema", "write", "--for", "all", "--output", "schemas", "--dry-run"],
    );
    let expected = indoc! {"
      schemas/version-info.schema.json
      schemas/progress-event.schema.json
      schemas/error-response.schema.json
      schemas/input-config-pipeline.schema.json
      schemas/input-config-timetree.schema.json
      schemas/input-config-optimize.schema.json
      schemas/input-config-prune.schema.json
      schemas/input-config-ancestral.schema.json
      schemas/input-config-homoplasy.schema.json
      schemas/input-config-clock.schema.json
      schemas/input-config-mugration.schema.json
    "};
    assert_eq!(
      (
        Run {
          code: Some(0),
          stdout: expected.to_owned(),
          stderr: String::new(),
        },
        false
      ),
      (run, dir.path().join("schemas").exists())
    );
  }

  #[test]
  fn test_cli_schema_write_dry_run_of_one_schema_names_standard_output() {
    let dir = tempdir().unwrap();
    let run = treetime_in(dir.path(), &["schema", "write", "--for", "timetree", "--dry-run"]);
    assert_eq!(
      Run {
        code: Some(0),
        stdout: "-\n".to_owned(),
        stderr: String::new(),
      },
      run
    );
  }

  mod helpers {
    use std::fs;
    use std::path::Path;
    use std::process::Command;

    pub(crate) const TREE: &str = "((A:0.01,B:0.02):0.01,C:0.03);\n";

    pub(crate) const ALIGNMENT: &str = ">A\nACGTACGTAC\n>B\nACGTACGTAA\n>C\nACGAACGTAC\n";

    #[derive(Debug, PartialEq, Eq)]
    pub(crate) struct Run {
      pub(crate) code: Option<i32>,
      pub(crate) stdout: String,
      pub(crate) stderr: String,
    }

    pub(crate) fn treetime_in(dir: &Path, args: &[&str]) -> Run {
      let output = Command::new(env!("CARGO_BIN_EXE_treetime"))
        .args(args)
        .current_dir(dir)
        .env("NO_COLOR", "1")
        .env_remove("RUST_BACKTRACE")
        .env_remove("RUST_LIB_BACKTRACE")
        .output()
        .unwrap();
      Run {
        code: output.status.code(),
        stdout: String::from_utf8(output.stdout).unwrap(),
        stderr: String::from_utf8(output.stderr).unwrap(),
      }
    }

    pub(crate) fn planned_lines(stdout: &str) -> Vec<String> {
      stdout
        .lines()
        .skip_while(|line| *line != "# Planned output files:")
        .map(str::to_owned)
        .collect()
    }

    pub(crate) fn yaml_block(stdout: &str, key: &str) -> String {
      let header = format!("{key}:");
      stdout
        .lines()
        .skip_while(|line| *line != header)
        .enumerate()
        .take_while(|(index, line)| *index == 0 || line.starts_with([' ', '-']))
        .fold(String::new(), |mut block, (_, line)| {
          block.push_str(line);
          block.push('\n');
          block
        })
    }

    pub(crate) fn path_in(dir: &Path, relative: &str) -> String {
      dir.join(relative).to_str().unwrap().to_owned()
    }

    pub(crate) fn write_file(dir: &Path, name: &str, text: &str) {
      fs::write(dir.join(name), text).unwrap();
    }
  }
}
