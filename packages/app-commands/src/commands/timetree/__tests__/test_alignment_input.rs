#[cfg(test)]
mod tests {
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use crate::commands::timetree::run::run_timetree_estimation;
  use clap::Parser;
  use eyre::Report;
  use std::fs;
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_utils::assert_error;

  #[test]
  fn test_timetree_one_third_rule_offers_no_override() -> Result<(), Report> {
    let dir = tempdir()?;
    let tree = dir.path().join("tree.nwk");
    let alignment = dir.path().join("aln.fasta");
    let metadata = dir.path().join("metadata.tsv");
    fs::write(&tree, "(A:0.01,B:0.02,C:0.03)root;\n")?;
    fs::write(&alignment, ">A\nACGT\n")?;
    fs::write(&metadata, "name\tdate\nA\t2000\nB\t2001\nC\t2002\n")?;
    let raw = TreetimeTimetreeArgsRaw::try_parse_from([
      "timetree",
      &format!("--tree={}", tree.display()),
      &format!("--alignment={}", alignment.display()),
      &format!("--metadata={}", metadata.display()),
      &format!("--output-all={}", dir.path().display()),
      "--seed=7",
    ])?;
    let args = TreetimeTimetreeArgs::try_from(raw)?;

    assert_error!(
      run_timetree_estimation(&args, &NoopCancel, &NoopProgress, &NoopProgress),
      format!(
        "In alignment '{}': More than one third of terminal nodes (2 of 3) cannot be assigned a sequence. \
         Are you sure the alignment belongs to the tree?",
        alignment.display()
      )
    );
    Ok(())
  }
}
