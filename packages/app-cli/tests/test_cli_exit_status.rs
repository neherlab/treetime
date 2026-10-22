#[cfg(test)]
mod tests {
  use helpers::{ESC, Run, data_path, treetime, write_file};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use tempfile::tempdir;

  #[test]
  fn test_cli_missing_input_is_a_plain_usage_error() {
    let run = treetime(&["ancestral"], &[]);
    assert_eq!(
      Run {
        code: Some(2),
        stdout: String::new(),
        stderr: "error: the following required arguments were not provided:\n  --tree <TREE>\n".to_owned(),
      },
      run
    );
  }

  #[test]
  fn test_cli_argument_check_after_parsing_is_a_plain_usage_error() {
    let run = treetime(
      &["timetree", "--tree", "tree.nwk", "--branch-length-mode", "marginal"],
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
  fn test_cli_run_without_outputs_is_a_plain_usage_error() {
    let run = treetime(
      &[
        "clock",
        "--tree",
        &data_path("zika/20/tree.nwk"),
        "--metadata",
        &data_path("zika/20/metadata.tsv"),
      ],
      &[],
    );
    assert_eq!(
      Run {
        code: Some(2),
        stdout: String::new(),
        stderr: "error: No output flags provided. At least one is required: --output-all or one of the \
                 --output-tree-* / --output-* flags\n"
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
  fn test_cli_failed_run_exits_with_one_and_shows_no_source_location() {
    let run = treetime(
      &[
        "clock",
        "--tree",
        "missing/tree.nwk",
        "--metadata",
        "missing/metadata.tsv",
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
  #[case::failure_never(   &["clock", "--tree=x.nwk", "--metadata=x.tsv", "--color=never"],  &[],        (Some(1), false))]
  #[case::failure_always(  &["clock", "--tree=x.nwk", "--metadata=x.tsv", "--color=always"], &[],        (Some(1), true))]
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
