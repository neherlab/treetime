#[cfg(test)]
mod tests {
  use crate::commands::prune::args::{TreetimePruneArgs, TreetimePruneArgsRaw};
  use crate::commands::prune::run::run_prune;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::output_args::{NewickAnnotationsArg, OutputCoreArgs};
  use crate::commands::shared::tree_input::TreeArgs;
  use eyre::{Report, WrapErr};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_utils::io::fs::read_file_to_string;

  const TREE: &str = "((A:0.1,B:0.1)AB:0.1,(C:0.1,D:0.1)CD:0.1)root;";

  #[rustfmt::skip]
  #[rstest]
  #[case::substitution( ">A\nACGT\n>B\nACGT\n>C\nTCGT\n>D\nACGT\n", "((A:0.1,B:0.1)AB:0.1,(C[&mutations=\"A1T\"]:0.1,D:0.1)CD:0.1)root;")]
  #[case::terminal_gap( ">A\nACGT\n>B\nACGT\n>C\n-CGT\n>D\nACGT\n", "((A:0.1,B:0.1)AB:0.1,(C:0.1,D:0.1)CD:0.1)root;")]
  #[trace]
  fn test_prune_tree_carries_mutations_of_the_alignment_without_a_pruning_operation(
    #[case] alignment: &str,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let dir = tempdir().wrap_err("When creating a temporary directory")?;
    let tree_path = dir.path().join("tree.nwk");
    let fasta_path = dir.path().join("aln.fasta");
    let nwk_path = dir.path().join("out.nwk");
    std::fs::write(&tree_path, TREE).wrap_err("When writing the tree fixture")?;
    std::fs::write(&fasta_path, alignment).wrap_err("When writing the alignment fixture")?;

    let args = TreetimePruneArgs::try_from(TreetimePruneArgsRaw {
      alignment: AlignmentArgs {
        alignment: vec![fasta_path],
      },
      tree: TreeArgs { tree: Some(tree_path) },
      output: OutputCoreArgs {
        newick_annotations: NewickAnnotationsArg::Beast,
        output_tree_nwk: Some(nwk_path.clone()),
        ..OutputCoreArgs::default()
      },
      ..TreetimePruneArgsRaw::default()
    })?;
    run_prune(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

    assert_eq!(expected, read_file_to_string(&nwk_path)?.trim());
    Ok(())
  }
}
