#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{RecordingLog, project_root};
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::tree_input::TreeArgs;
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use crate::commands::timetree::branch_length_mode::BranchLengthModeCli;
  use crate::commands::timetree::initialization::load_input_data;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::path::PathBuf;
  use treetime::progress::LogLevel;
  use treetime::timetree::params::BranchLengthMode;
  use treetime_utils::vec_of_owned;

  #[rustfmt::skip]
  #[rstest]
  #[case::input(                  (BranchLengthModeCli::Input,    true,  Some(0.01)), BranchLengthMode::Input)]
  #[case::marginal(               (BranchLengthModeCli::Marginal, true,  Some(0.5)),  BranchLengthMode::Marginal)]
  #[case::auto_without_alignment( (BranchLengthModeCli::Auto,     false, Some(0.01)), BranchLengthMode::Input)]
  #[case::auto_short_branches(    (BranchLengthModeCli::Auto,     true,  Some(0.05)), BranchLengthMode::Marginal)]
  #[case::auto_at_the_limit(      (BranchLengthModeCli::Auto,     true,  Some(0.1)),  BranchLengthMode::Marginal)]
  #[case::auto_long_branches(     (BranchLengthModeCli::Auto,     true,  Some(0.15)), BranchLengthMode::Input)]
  #[case::auto_without_lengths(   (BranchLengthModeCli::Auto,     true,  None),       BranchLengthMode::Marginal)]
  #[trace]
  fn test_branch_length_mode_resolves(
    #[case] (mode, has_alignment, longest_branch): (BranchLengthModeCli, bool, Option<f64>),
    #[case] expected: BranchLengthMode,
  ) {
    assert_eq!(expected, mode.resolve(has_alignment, longest_branch));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::auto_without_alignment( (BranchLengthModeCli::Auto,     "flu/h3n2/20", false), (BranchLengthMode::Input,    vec_of_owned!["Branch length mode auto: no alignment, using input"]))]
  #[case::auto_short_branches(    (BranchLengthModeCli::Auto,     "flu/h3n2/20", true),  (BranchLengthMode::Marginal, vec_of_owned!["Branch length mode auto: the longest input branch is 2.186e-2 substitutions per site, using marginal"]))]
  #[case::auto_long_branches(     (BranchLengthModeCli::Auto,     "tb/20",       true),  (BranchLengthMode::Input,    vec_of_owned!["Branch length mode auto: the longest input branch is 1.509e-1 substitutions per site, using input"]))]
  #[case::marginal_long_branches( (BranchLengthModeCli::Marginal, "tb/20",       true),  (BranchLengthMode::Marginal, vec![]))]
  #[trace]
  fn test_branch_length_mode_is_resolved_from_the_input_tree(
    #[case] (mode, dataset, with_alignment): (BranchLengthModeCli, &str, bool),
    #[case] (expected_mode, expected_messages): (BranchLengthMode, Vec<String>),
  ) -> Result<(), Report> {
    let dir = project_root().join("data").join(dataset);
    let alignment = with_alignment.then(|| dir.join("aln.fasta.xz"));

    let (mode, messages) = helpers::resolved_mode(mode, dir.join("tree.nwk"), alignment)?;

    assert_eq!((expected_mode, expected_messages), (mode, messages));
    Ok(())
  }

  #[test]
  fn test_branch_length_mode_auto_without_branch_lengths_is_marginal() -> Result<(), Report> {
    let dir = tempfile::tempdir()?;
    let tree = dir.path().join("tree.nwk");
    let alignment = dir.path().join("aln.fasta");
    std::fs::write(&tree, "((A,B),C);")?;
    std::fs::write(&alignment, ">A\nACGT\n>B\nACGA\n>C\nTCGA\n")?;

    let (mode, messages) = helpers::resolved_mode(BranchLengthModeCli::Auto, tree, Some(alignment))?;

    assert_eq!(
      (
        BranchLengthMode::Marginal,
        vec_of_owned!["Branch length mode auto: no input branch has a length, using marginal"]
      ),
      (mode, messages)
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn resolved_mode(
      branch_length_mode: BranchLengthModeCli,
      tree: PathBuf,
      alignment: Option<PathBuf>,
    ) -> Result<(BranchLengthMode, Vec<String>), Report> {
      let sequence_length = alignment.is_none().then_some(1000);
      let args = TreetimeTimetreeArgs::try_from(TreetimeTimetreeArgsRaw {
        alignment: AlignmentArgs {
          alignment: alignment.into_iter().collect(),
        },
        tree: TreeArgs { tree: Some(tree) },
        sequence_length,
        branch_length_mode,
        ..TreetimeTimetreeArgsRaw::default()
      })?;
      let log = RecordingLog::default();

      let input = load_input_data(&args, &log)?;

      let messages = log
        .messages(LogLevel::Info)
        .into_iter()
        .filter(|message| message.starts_with("Branch length mode"))
        .collect();
      Ok((input.branch_length_mode, messages))
    }
  }
}
