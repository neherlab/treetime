#[cfg(test)]
mod tests {
  use helpers::{Run, data_path, treetime};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use tempfile::tempdir;

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree( &["--tree", "<tree>", "--aln", "<aln>", "--dates", "<metadata>", "--outdir", "<out>", "--rng-seed", "1"],      "timetree.nwk")]
  #[case::ancestral(&["ancestral", "--aln", "<aln>", "--tree", "<tree>", "--gtr", "JC69", "--marginal", "--outdir", "<out>"],      "ancestral.nwk")]
  #[case::clock(    &["clock", "--tree", "<tree>", "--dates", "<metadata>", "--aln", "<aln>", "--reroot", "best", "--outdir", "<out>"], "clock.nwk")]
  #[case::homoplasy(&["homoplasy", "--aln", "<aln>", "--tree", "<tree>", "--gtr", "JC69", "--outdir", "<out>"],                    "homoplasy.nwk")]
  #[case::mugration(&["mugration", "--tree", "<tree>", "--states", "<metadata>", "--attribute", "country", "--confidence", "--outdir", "<out>"],
                                                                                                                                        "mugration.nwk")]
  #[trace]
  fn test_cli_v0_spellings_run_a_typical_v0_command_line(#[case] args: &[&str], #[case] output: &str) {
    let dir = tempdir().unwrap();
    let out = dir.path().join("out");
    let args: Vec<String> = args
      .iter()
      .map(|arg| {
        arg
          .replace("<tree>", &data_path("zika/20/tree.nwk"))
          .replace("<aln>", &data_path("zika/20/aln.fasta.xz"))
          .replace("<metadata>", &data_path("zika/20/metadata.tsv"))
          .replace("<out>", out.to_str().unwrap())
      })
      .collect();
    let args: Vec<&str> = args.iter().map(String::as_str).collect();

    let run = treetime(&args);

    assert_eq!((Some(0), true), (run.code, out.join(output).is_file()), "{}", run.stderr);
  }

  #[test]
  fn test_cli_v0_spellings_log_the_rewrites_at_debug_level() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("out");
    let tree = data_path("zika/20/tree.nwk");
    let aln = data_path("zika/20/aln.fasta.xz");
    let aln_arg = format!("--aln {aln}");

    let run = treetime(&[
      "ancestral",
      "--tree",
      &tree,
      "--aln",
      &aln,
      "--outdir",
      out.to_str().unwrap(),
      "--verbose",
      "4",
      "--no-progress",
    ]);

    let message = format!("TreeTime 0.x argument '{aln_arg}' read as '--alignment {aln}'");
    assert_eq!(
      (Some(0), true),
      (run.code, run.stderr.contains(&message)),
      "{}",
      run.stderr
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::removed_flag(  &["--tree", "t.nwk", "--covariation"],
                         "error: --covariation is not available in TreeTime 1: TreeTime 1 fits the root-to-tip regression by ordinary least squares, without the covariation model\n")]
  #[case::changed_value( &["timetree", "--time-marginal", "always"],
                         "error: --time-marginal=always is not available in TreeTime 1: use --confidence\n")]
  #[case::arg_command(   &["arg"],
                         "error: the TreeTime 0.x command 'arg' is not available in TreeTime 1\n")]
  #[trace]
  fn test_cli_v0_spellings_stop_with_a_usage_error(#[case] args: &[&str], #[case] stderr: &str) {
    let run = treetime(args);
    assert_eq!(
      Run {
        code: Some(2),
        stdout: String::new(),
        stderr: stderr.to_owned(),
      },
      run
    );
  }

  #[test]
  fn test_cli_v0_spellings_version_command_prints_the_version() {
    assert_eq!(treetime(&["--version"]), treetime(&["version"]));
  }

  mod helpers {
    use std::path::Path;
    use std::process::Command;

    #[derive(Debug, PartialEq, Eq)]
    pub(crate) struct Run {
      pub(crate) code: Option<i32>,
      pub(crate) stdout: String,
      pub(crate) stderr: String,
    }

    pub(crate) fn treetime(args: &[&str]) -> Run {
      let output = Command::new(env!("CARGO_BIN_EXE_treetime"))
        .args(args)
        .env_remove("NO_COLOR")
        .env_remove("CLICOLOR")
        .env_remove("CLICOLOR_FORCE")
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

    pub(crate) fn data_path(relative: &str) -> String {
      Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data")
        .join(relative)
        .to_str()
        .unwrap()
        .to_owned()
    }
  }
}
