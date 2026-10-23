#[cfg(test)]
mod tests {
  use helpers::{ESC, Run, data_path, treetime, write_file};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::fs;
  use std::path::Path;
  use tempfile::tempdir;
  use treetime_utils::o;

  #[test]
  fn test_cli_missing_input_is_a_plain_usage_error() {
    let run = treetime(&["ancestral"], &[]);
    assert_eq!(
      Run {
        code: Some(2),
        stdout: String::new(),
        stderr: "error: the following required arguments were not provided:\n  --tree <FILE>\n".to_owned(),
      },
      run
    );
  }

  #[test]
  fn test_cli_argument_check_after_parsing_is_a_plain_usage_error() {
    let run = treetime(
      &[
        "timetree",
        "--tree",
        "tree.nwk",
        "--metadata",
        "metadata.tsv",
        "--branch-length-mode",
        "marginal",
      ],
      &[],
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
  fn test_cli_timetree_without_metadata_is_a_plain_usage_error() {
    let run = treetime(&["timetree", "--tree", "tree.nwk", "--sequence-length", "1000"], &[]);
    assert_eq!(
      Run {
        code: Some(2),
        stdout: String::new(),
        stderr: "error: the following required arguments were not provided:\n  --metadata <FILE>\n".to_owned(),
      },
      run
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::hky85_aa(&["--alphabet", "aa", "--model", "hky85"], "the model hky85 needs the alphabet nuc because it sets rates between nucleotides, but the alphabet is aa")]
  #[case::jtt92_aa(&["--alphabet", "aa", "--model", "jtt92"], "the model jtt92 needs the alphabet aa-no-stop because its matrix holds the 20 amino acids without the stop '*', but the alphabet is aa")]
  #[trace]
  fn test_cli_timetree_without_alignment_checks_the_model_against_the_alphabet(
    #[case] model: &[&str],
    #[case] message: &str,
  ) {
    let dir = tempdir().unwrap();
    let (tree, metadata) = (data_path("zika/20/tree.nwk"), data_path("zika/20/metadata.tsv"));
    let output = dir.path().to_str().unwrap();
    let args = [
      &["timetree", "--tree", &tree, "--metadata", &metadata, "--sequence-length", "1000", "--output-all", output][..],
      model,
    ]
    .concat();
    assert_eq!(
      Run {
        code: Some(2),
        stdout: String::new(),
        stderr: format!("error: {message}\n"),
      },
      treetime(&args, &[])
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::ancestral(&["ancestral", "--tree", "missing.nwk", "--alignment", "missing.fasta"])]
  #[case::clock(    &["clock",     "--tree", "missing.nwk", "--metadata", "missing.tsv"])]
  #[case::homoplasy(&["homoplasy", "--tree", "missing.nwk", "--alignment", "missing.fasta"])]
  #[case::mugration(&["mugration", "--tree", "missing.nwk", "--metadata", "missing.tsv", "--attribute", "country"])]
  #[case::optimize( &["optimize",  "--tree", "missing.nwk", "--alignment", "missing.fasta"])]
  #[case::prune(    &["prune",     "--tree", "missing.nwk"])]
  #[case::timetree( &["timetree",  "--tree", "missing.nwk", "--metadata", "missing.tsv", "--alignment", "missing.fasta"])]
  #[trace]
  fn test_cli_run_without_outputs_fails_before_reading_inputs(#[case] args: &[&str]) {
    let run = treetime(args, &[]);
    assert_eq!(
      Run {
        code: Some(2),
        stdout: String::new(),
        stderr: "error: No outputs selected: add -O <DIR> to write the default outputs into a folder, or a per-file \
                 --output-* flag\n"
          .to_owned(),
      },
      run
    );
  }

  #[test]
  fn test_cli_output_selection_without_output_folder_is_a_plain_usage_error() {
    let dir = tempdir().unwrap();
    let config = write_file(dir.path(), "config.yaml", "output_selection: [nwk]\n");
    let run = treetime(
      &[
        "clock",
        "--config",
        &config,
        "--tree",
        &data_path("zika/20/tree.nwk"),
        "--metadata",
        &data_path("zika/20/metadata.tsv"),
        "--output-tree-nexus",
        dir.path().join("tree.nexus").to_str().unwrap(),
      ],
      &[],
    );
    assert_eq!(
      Run {
        code: Some(2),
        stdout: String::new(),
        stderr: "error: --output-selection requires --output-all\n".to_owned(),
      },
      run
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::conflict(    "clock", "keep_root: true\nreroot: oldest\n",   "error: invalid configuration: `keep_root` cannot be used together with `reroot`")]
  #[case::unknown_key( "clock", "definitely_not_a_real_field: 1\n",   "error: invalid configuration: unknown field `definitely_not_a_real_field`")]
  #[case::partition_and_model("ancestral", "partition:\n  - {name: all}\nmodel: jc69\n", "error: invalid configuration: `model` cannot be used together with `partition`")]
  #[trace]
  fn test_cli_invalid_config_file_is_a_usage_error(#[case] command: &str, #[case] text: &str, #[case] first_line: &str) {
    let dir = tempdir().unwrap();
    let config = write_file(dir.path(), "config.yaml", text);
    let run = treetime(&[command, "--config", &config], &[]);
    assert_eq!(
      (Some(2), first_line, false, false),
      (
        run.code,
        run.stderr.lines().next().unwrap_or_default(),
        run.stderr.contains(ESC),
        run.stderr.contains("Location:")
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::command(  "timetree", "coalescent: 1e999y\n",                                      "")]
  #[case::pipeline( "pipeline", "steps:\n  - name: tt\n    timetree:\n      coalescent: 1e999y\n", "in pipeline step `tt`: ")]
  #[trace]
  fn test_cli_config_value_that_passes_the_schema_and_fails_conversion_is_a_usage_error(
    #[case] command: &str,
    #[case] text: &str,
    #[case] context: &str,
  ) {
    let dir = tempdir().unwrap();
    let config = write_file(dir.path(), "config.yaml", text);
    let run = treetime(&[command, "--config", &config], &[]);
    let first_line = format!(
      "error: invalid configuration: {context}invalid value: '1e999y' is not a coalescent prior. Give a fixed time \
       scale as a positive number in substitutions per site (0.01, 1e-2) or as a positive number directly followed by \
       the time unit y, w or d (10y, 26w, 90d); opt to fit a constant time scale; or skyline to fit a \
       piecewise-constant one"
    );
    assert_eq!(
      (Some(2), first_line.as_str(), false),
      (
        run.code,
        run.stderr.lines().next().unwrap_or_default(),
        run.stderr.contains("Location:")
      )
    );
  }

  #[test]
  fn test_cli_run_without_outputs_fails_before_reading_a_malformed_input() {
    let dir = tempdir().unwrap();
    let tree = write_file(dir.path(), "tree.nwk", "this is not a Newick tree");
    let run = treetime(&["prune", "--tree", &tree], &[]);
    assert_eq!(
      Run {
        code: Some(2),
        stdout: String::new(),
        stderr: "error: No outputs selected: add -O <DIR> to write the default outputs into a folder, or a per-file \
                 --output-* flag\n"
          .to_owned(),
      },
      run
    );
  }

  #[test]
  fn test_cli_run_that_would_overwrite_its_input_fails_before_writing() {
    let dir = tempdir().unwrap();
    let tree = write_file(dir.path(), "tree.nwk", "(A:0.1,B:0.2)root;\n");
    let run = treetime(&["prune", "--tree", &tree, "--output-tree-nwk", &tree], &[]);
    assert_eq!(
      (
        Run {
          code: Some(2),
          stdout: String::new(),
          stderr: format!(
            "error: --tree '{tree}' names the same file as the output --output-tree-nwk '{tree}': the run would \
             overwrite its own input\n"
          ),
        },
        o!("(A:0.1,B:0.2)root;\n")
      ),
      (run, fs::read_to_string(&tree).unwrap())
    );
  }

  #[test]
  fn test_cli_translation_template_that_matches_an_output_fails_before_writing() {
    let dir = tempdir().unwrap();
    let tree = write_file(dir.path(), "tree.nwk", "(A:0.1,B:0.2)root;\n");
    let alignment = write_file(dir.path(), "aln.fasta", ">A\nACGT\n>B\nACGA\n");
    let template = dir.path().join("aa/{cds}.fasta").to_str().unwrap().to_owned();
    let output = dir.path().join("aa/S.fasta").to_str().unwrap().to_owned();
    let run = treetime(
      &[
        "ancestral",
        "--tree",
        &tree,
        "--alignment",
        &alignment,
        "--translations",
        &template,
        "--cds-selection",
        "S",
        "--output-sequences",
        &output,
      ],
      &[],
    );
    assert_eq!(
      (
        Run {
          code: Some(2),
          stdout: String::new(),
          stderr: format!(
            "error: --translations '{template}' can name the same file as the output --output-sequences '{output}': \
             the run would overwrite its own input\n"
          ),
        },
        false
      ),
      (run, Path::new(&output).exists())
    );
  }

  #[test]
  fn test_cli_homoplasy_without_alignment_fails_before_reading_the_drms() {
    let dir = tempdir().unwrap();
    let tree = write_file(dir.path(), "tree.nwk", "(A:0.1,B:0.2)root;\n");
    let drms = dir.path().join("missing.tsv");
    let out = dir.path().join("out");
    let run = treetime(
      &[
        "homoplasy",
        "--tree",
        &tree,
        "--drms",
        drms.to_str().unwrap(),
        "-O",
        out.to_str().unwrap(),
      ],
      &[],
    );
    assert_eq!(
      Run {
        code: Some(2),
        stdout: String::new(),
        stderr: "error: --alignment is required: pass one or more FASTA files, or '-' to read standard input\n"
          .to_owned(),
      },
      run
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::generated_output_overwrites_input(
    "output_all: out\nsteps:\n  - name: anc\n    ancestral:\n      alignment: [a.fasta]\n      tree: out/anc/ancestral.nwk\n",
    "error: in pipeline step `anc`: --tree 'DIR/out/anc/ancestral.nwk' names the same file as the output --output-all \
     'DIR/out/anc/ancestral.nwk': the run would overwrite its own input\n"
  )]
  #[case::step_without_outputs(
    "steps:\n  - name: anc\n    ancestral:\n      alignment: [a.fasta]\n      tree: out/anc/ancestral.nwk\n",
    "error: in pipeline step `anc`: the step writes no outputs; set the top-level `output_all` key of the pipeline, \
     pass `treetime pipeline --output-all <dir>`, or name output files in the step\n"
  )]
  #[trace]
  fn test_cli_pipeline_step_output_error_is_a_usage_error(#[case] pipeline: &str, #[case] expected_stderr: &str) {
    let dir = tempdir().unwrap();
    fs::create_dir_all(dir.path().join("out/anc")).unwrap();
    let tree = write_file(&dir.path().join("out/anc"), "ancestral.nwk", "(A:0.1,B:0.2)root;\n");
    let config = write_file(dir.path(), "pipeline.yaml", pipeline);
    let run = treetime(&["pipeline", "--config", &config], &[]);
    assert_eq!(
      (
        Run {
          code: Some(2),
          stdout: String::new(),
          stderr: expected_stderr.replace("DIR", dir.path().to_str().unwrap()),
        },
        o!("(A:0.1,B:0.2)root;\n")
      ),
      (run, fs::read_to_string(&tree).unwrap())
    );
  }

  #[test]
  fn test_cli_standard_input_for_two_inputs_is_a_usage_error() {
    let dir = tempdir().unwrap();
    let output = dir.path().join("out.nwk");
    let run = treetime(
      &[
        "prune",
        "--tree",
        "-",
        "--remove-file",
        "-",
        "--output-tree-nwk",
        output.to_str().unwrap(),
      ],
      &[],
    );
    assert_eq!(
      Run {
        code: Some(2),
        stdout: String::new(),
        stderr: "error: Standard input can be read only once, but --tree and --remove-file read it\n".to_owned(),
      },
      run
    );
  }

  #[test]
  fn test_cli_failed_run_exits_with_one_and_shows_no_source_location() {
    let run = treetime(
      &[
        "clock",
        "--tree",
        "missing/tree.nwk",
        "--metadata",
        "missing/metadata.tsv",
        "--output-tree-nwk",
        "missing/tree.out.nwk",
      ],
      &[],
    );
    assert_eq!(
      (Some(1), true, false, false, false),
      (
        run.code,
        run.stderr.starts_with("Error:"),
        run.stderr.contains(ESC),
        run.stderr.contains("Location:"),
        run.stderr.contains("RUST_BACKTRACE")
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::usage_never(     &["ancestral", "--color=never"],                          &[],                (Some(2), false))]
  #[case::usage_always(    &["ancestral", "--color=always"],                         &[],                (Some(2), true))]
  #[case::usage_auto_pipe( &["ancestral"],                                           &[],                (Some(2), false))]
  #[case::no_color_always( &["ancestral", "--color", "always"],                      &[("NO_COLOR", "1")], (Some(2), true))]
  #[case::clap_never(      &["ancestral", "--no-such-flag", "--color=never"],        &[],                (Some(2), false))]
  #[case::clap_always(     &["--color=always", "ancestral", "--no-such-flag"],       &[],                (Some(2), true))]
  #[case::failure_never(   &["clock", "--tree=x.nwk", "--metadata=x.tsv", "--output-tree-nwk=x/out.nwk", "--color=never"],  &[], (Some(1), false))]
  #[case::failure_always(  &["clock", "--tree=x.nwk", "--metadata=x.tsv", "--output-tree-nwk=x/out.nwk", "--color=always"], &[], (Some(1), true))]
  #[trace]
  fn test_cli_color_flag_controls_error_colors(
    #[case] args: &[&str],
    #[case] env: &[(&str, &str)],
    #[case] expected: (Option<i32>, bool),
  ) {
    let run = treetime(args, env);
    assert_eq!(expected, (run.code, run.stderr.contains(ESC)));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::never(  "never",  false)]
  #[case::always( "always", true)]
  #[trace]
  fn test_cli_color_flag_controls_log_colors(#[case] color: &str, #[case] expected_colors: bool) {
    let dir = tempdir().unwrap();
    let out = dir.path().join("out");
    let run = treetime(
      &[
        "clock",
        "--tree",
        &data_path("zika/20/tree.nwk"),
        "--metadata",
        &data_path("zika/20/metadata.tsv"),
        "--output-all",
        out.to_str().unwrap(),
        "--verbose",
        "--no-progress",
        "--color",
        color,
      ],
      &[],
    );
    assert_eq!(
      (Some(0), true, expected_colors),
      (run.code, run.stderr.contains("# Command line arguments"), run.stderr.contains(ESC))
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::analysis_command( &["clock", "--jobs", "1", "--help"],  Some(0))]
  #[case::pipeline(         &["pipeline", "-j", "1", "--help"],   Some(0))]
  #[case::before_command(   &["--jobs", "1", "clock", "--help"],  Some(2))]
  #[case::tool_command(     &["completions", "bash", "-j", "1"],  Some(2))]
  #[trace]
  fn test_cli_jobs_flag_belongs_to_the_run_commands(#[case] args: &[&str], #[case] expected: Option<i32>) {
    assert_eq!(expected, treetime(args, &[]).code);
  }

  #[test]
  fn test_cli_schema_write_prints_the_config_schema() {
    let run = treetime(&["schema", "write", "--for", "timetree", "-o", "-"], &[]);
    let schema: serde_json::Value = serde_json::from_str(&run.stdout).unwrap();
    assert_eq!((Some(0), Some("TimetreeConfig")), (run.code, schema["title"].as_str()));
  }

  #[test]
  fn test_cli_pipeline_dry_run_prints_the_plan_and_writes_nothing() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("out");
    let run = treetime(
      &[
        "pipeline",
        "--config",
        &data_path("zika/20/pipeline.yaml"),
        "--dry-run",
        "-O",
        out.to_str().unwrap(),
      ],
      &[],
    );
    assert_eq!(
      (Some(0), true, false),
      (run.code, run.stdout.contains(out.to_str().unwrap()), out.exists())
    );
  }

  mod helpers {
    use std::fs;
    use std::path::Path;
    use std::process::Command;

    pub(crate) const ESC: char = '\u{1b}';

    #[derive(Debug, PartialEq, Eq)]
    pub(crate) struct Run {
      pub(crate) code: Option<i32>,
      pub(crate) stdout: String,
      pub(crate) stderr: String,
    }

    pub(crate) fn treetime(args: &[&str], env: &[(&str, &str)]) -> Run {
      let output = Command::new(env!("CARGO_BIN_EXE_treetime"))
        .args(args)
        .env_remove("NO_COLOR")
        .env_remove("CLICOLOR")
        .env_remove("CLICOLOR_FORCE")
        .env_remove("RUST_BACKTRACE")
        .env_remove("RUST_LIB_BACKTRACE")
        .envs(env.iter().copied())
        .output()
        .unwrap();
      Run {
        code: output.status.code(),
        stdout: String::from_utf8(output.stdout).unwrap(),
        stderr: String::from_utf8(output.stderr).unwrap(),
      }
    }

    pub(crate) fn data_path(relative: &str) -> String {
      Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data")
        .join(relative)
        .to_str()
        .unwrap()
        .to_owned()
    }

    pub(crate) fn write_file(dir: &Path, name: &str, text: &str) -> String {
      let path = dir.join(name);
      fs::write(&path, text).unwrap();
      path.to_str().unwrap().to_owned()
    }
  }
}
