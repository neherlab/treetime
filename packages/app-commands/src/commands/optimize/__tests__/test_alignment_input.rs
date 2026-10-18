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
  use treetime_utils::{assert_error, o};

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

    helpers::run(dir.path())?;

    let tree = nwk_read_file(dir.path().join("optimize.nwk"))?;
    let names = tree.names();
    let leaves = tree
      .graph
      .get_leaves()
      .map(|leaf| names[&leaf.key()].clone())
      .sorted()
      .collect_vec();
    assert_eq!(vec![Some(o!("A")), Some(o!("B")), Some(o!("C")), Some(o!("D"))], leaves);
    Ok(())
  }

  #[test]
  fn test_optimize_one_third_rule_offers_no_override() -> Result<(), Report> {
    let dir = tempdir()?;
    fs::write(dir.path().join("tree.nwk"), "(A:0.01,B:0.02,C:0.03)root;\n")?;
    fs::write(dir.path().join("aln.fasta"), ">A\nACGT\n")?;

    assert_error!(
      helpers::run(dir.path()),
      format!(
        "In alignment '{}': More than one third of terminal nodes (2 of 3) cannot be assigned a sequence. \
         Are you sure the alignment belongs to the tree?",
        dir.path().join("aln.fasta").display()
      )
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn run(dir: &Path) -> Result<(), Report> {
      let raw = TreetimeOptimizeArgsRaw::try_parse_from([
        "optimize",
        &format!("--tree={}", dir.join("tree.nwk").display()),
        &format!("--alignment={}", dir.join("aln.fasta").display()),
        &format!("--output-all={}", dir.display()),
        "--max-iter=1",
      ])?;
      let args = TreetimeOptimizeArgs::try_from(raw)?;
      run_optimize(&args, &NoopCancel, &NoopProgress, &NoopProgress)
    }
  }
}
