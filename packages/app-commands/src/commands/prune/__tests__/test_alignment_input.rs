#[cfg(test)]
mod tests {
  use crate::commands::prune::args::{TreetimePruneArgs, TreetimePruneArgsRaw};
  use crate::commands::prune::run::run_prune;
  use clap::Parser;
  use eyre::Report;
  use std::fs;
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_utils::assert_error;

  #[test]
  fn test_prune_one_third_rule_offers_no_override() -> Result<(), Report> {
    let dir = tempdir()?;
    let tree = dir.path().join("tree.nwk");
    let alignment = dir.path().join("aln.fasta");
    fs::write(&tree, "(A:0.01,B:0.02,C:0.03)root;\n")?;
    fs::write(&alignment, ">A\nACGT\n")?;
    let raw = TreetimePruneArgsRaw::try_parse_from([
      "prune",
      &format!("--tree={}", tree.display()),
      &format!("--alignment={}", alignment.display()),
      &format!("--output-all={}", dir.path().display()),
      "--prune-empty",
    ])?;
    let args = TreetimePruneArgs::try_from(raw)?;

    assert_error!(
      run_prune(&args, &NoopCancel, &NoopProgress, &NoopProgress),
      format!(
        "In alignment '{}': More than one third of terminal nodes (2 of 3) cannot be assigned a sequence. \
         Are you sure the alignment belongs to the tree?",
        alignment.display()
      )
    );
    Ok(())
  }
}
