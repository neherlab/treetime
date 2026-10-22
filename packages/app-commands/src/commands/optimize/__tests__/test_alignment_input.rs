#[cfg(test)]
mod tests {
  use crate::commands::optimize::args::{TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
  use crate::commands::optimize::run::run_optimize;
  use clap::Parser;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::fs;
  use std::path::Path;
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_io::nwk::nwk_read_file;
  use treetime_utils::{assert_error, o, vec_of_owned};

  #[test]
  fn test_optimize_fills_a_leaf_without_a_sequence() -> Result<(), Report> {
    let dir = tempdir()?;
    fs::write(
      dir.path().join("tree.nwk"),
      "((A:0.01,B:0.02):0.01,(C:0.01,D:0.02):0.01)root;\n",
    )?;
    fs::write(
      dir.path().join("aln.fasta"),
      ">A\nACGTACGTAC\n>B\nACGTACGTTC\n>C\nACGAACGTAC\n",
    )?;

    helpers::run(dir.path(), &[])?;

    let tree = nwk_read_file(dir.path().join("optimize.nwk"))?;
    let leaves = tree
      .graph
      .get_leaves()
      .map(|leaf| tree.names[leaf.key()].to_owned())
      .sorted()
      .collect_vec();
    assert_eq!(vec_of_owned!["A", "B", "C", "D"], leaves);
    Ok(())
  }

  #[test]
  fn test_optimize_one_third_rule_names_the_override() -> Result<(), Report> {
    let dir = tempdir()?;
    fs::write(dir.path().join("tree.nwk"), "(A:0.01,B:0.02,C:0.03)root;\n")?;
    fs::write(dir.path().join("aln.fasta"), ">A\nACGT\n")?;

    assert_error!(
      helpers::run(dir.path(), &[]),
      format!(
        "In alignment '{}': More than one third of terminal nodes (2 of 3) cannot be assigned a sequence. \
         Are you sure the alignment belongs to the tree? \
         Pass --allow-missing-sequences to proceed, treating missing tips as fully ambiguous.",
        dir.path().join("aln.fasta").display()
      )
    );
    Ok(())
  }

  #[test]
  fn test_optimize_allow_missing_sequences_lifts_the_one_third_rule() -> Result<(), Report> {
    let dir = tempdir()?;
    fs::write(dir.path().join("tree.nwk"), "(A:0.01,B:0.02,C:0.03)root;\n")?;
    fs::write(dir.path().join("aln.fasta"), ">A\nACGT\n")?;

    helpers::run(dir.path(), &["--allow-missing-sequences"])?;

    assert!(dir.path().join("optimize.nwk").exists());
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn run(dir: &Path, extra: &[&str]) -> Result<(), Report> {
      let args = [
        o!("optimize"),
        format!("--tree={}", dir.join("tree.nwk").display()),
        format!("--alignment={}", dir.join("aln.fasta").display()),
        format!("--output-all={}", dir.display()),
        o!("--max-iter=1"),
      ];
      let raw = TreetimeOptimizeArgsRaw::try_parse_from(args.iter().map(String::as_str).chain(extra.iter().copied()))?;
      let args = TreetimeOptimizeArgs::try_from(raw)?;
      run_optimize(&args, &NoopCancel, &NoopProgress, &NoopProgress)
    }
  }
}
