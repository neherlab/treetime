#[cfg(test)]
mod tests {
  use crate::commands::prune::args::{TreetimePruneArgs, TreetimePruneArgsRaw};
  use crate::commands::prune::run::run_prune;
  use clap::Parser;
  use eyre::Report;
  use indoc::indoc;
  use std::fs;
  use std::path::{Path, PathBuf};
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_prune_one_third_rule_names_the_override() -> Result<(), Report> {
    let dir = tempdir()?;
    let alignment = helpers::write_inputs(dir.path())?;

    assert_error!(
      helpers::run(dir.path(), &[]),
      format!(
        "In alignment '{}': More than one third of terminal nodes (2 of 3) cannot be assigned a sequence. \
         Are you sure the alignment belongs to the tree? \
         Pass --allow-missing-sequences to proceed, treating missing tips as fully ambiguous.",
        alignment.display()
      )
    );
    Ok(())
  }

  #[test]
  fn test_prune_allow_missing_sequences_lifts_the_one_third_rule() -> Result<(), Report> {
    let dir = tempdir()?;
    helpers::write_inputs(dir.path())?;

    helpers::run(dir.path(), &["--allow-missing-sequences"])
  }

  #[test]
  fn test_prune_amino_acid_alignment_runs() -> Result<(), Report> {
    let dir = tempdir()?;
    helpers::write_inputs(dir.path())?;
    fs::write(
      dir.path().join("aln.fasta"),
      indoc! {r"
        >A
        MKVL
        >B
        MKIL
        >C
        MRIL
      "},
    )?;

    helpers::run(dir.path(), &["--alphabet", "aa"])
  }

  mod helpers {
    use super::*;

    pub(super) fn write_inputs(dir: &Path) -> Result<PathBuf, Report> {
      let alignment = dir.join("aln.fasta");
      fs::write(dir.join("tree.nwk"), "(A:0.01,B:0.02,C:0.03)root;\n")?;
      fs::write(&alignment, ">A\nACGT\n")?;
      Ok(alignment)
    }

    pub(super) fn run(dir: &Path, extra: &[&str]) -> Result<(), Report> {
      let args = [
        o!("prune"),
        format!("--tree={}", dir.join("tree.nwk").display()),
        format!("--alignment={}", dir.join("aln.fasta").display()),
        format!("--output-all={}", dir.display()),
        o!("--merge-shared-mutations"),
      ];
      let raw = TreetimePruneArgsRaw::try_parse_from(args.iter().map(String::as_str).chain(extra.iter().copied()))?;
      let args = TreetimePruneArgs::try_from(raw)?;
      run_prune(&args, &NoopCancel, &NoopProgress, &NoopProgress)
    }
  }
}
