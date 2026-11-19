#[cfg(test)]
mod tests {
  use crate::commands::homoplasy::args::{TreetimeHomoplasyArgs, TreetimeHomoplasyArgsRaw};
  use clap::Parser;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::seq::mutation::PositionBase;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::defaults(&[],                                                                                      (0,   1.0, false, None,             10, false))]
  #[case::flags(   &["--const=100", "--rescale=0.5", "--detailed", "--drms=drms.tsv", "-n", "3", "--zero-based"], (100, 0.5, true,  Some("drms.tsv"), 3,  true))]
  #[trace]
  fn test_args_homoplasy_reads_homoplasy_flags(
    #[case] flags: &[&str],
    #[case] expected: (usize, f64, bool, Option<&str>, usize, bool),
  ) -> Result<(), Report> {
    let argv = [&"homoplasy", &"--tree=tree.nwk", &"--alignment=a.fasta"].into_iter().chain(flags);
    let args = TreetimeHomoplasyArgs::try_from(TreetimeHomoplasyArgsRaw::try_parse_from(argv)?)?;

    assert_eq!(
      expected,
      (
        args.constant_sites,
        args.rescale,
        args.detailed,
        args.drms.as_deref().and_then(|path| path.to_str()),
        args.num_mut,
        args.positions == PositionBase::Zero
      )
    );
    Ok(())
  }

  #[test]
  fn test_args_homoplasy_rejects_two_alignments() -> Result<(), Report> {
    let raw = TreetimeHomoplasyArgsRaw::try_parse_from([
      "homoplasy",
      "--tree=tree.nwk",
      "--alignment=a.fasta",
      "--alignment=b.fasta",
    ])?;
    assert_error!(
      TreetimeHomoplasyArgs::try_from(raw),
      "homoplasy reads exactly one --alignment file, but 2 were given: \
       the reported positions refer to the coordinates of one alignment"
    );
    Ok(())
  }

  #[test]
  fn test_args_homoplasy_requires_an_alignment() -> Result<(), Report> {
    let raw = TreetimeHomoplasyArgsRaw::try_parse_from(["homoplasy", "--tree=tree.nwk"])?;
    assert_error!(
      TreetimeHomoplasyArgs::try_from(raw),
      "--alignment is required: pass one or more FASTA files, or '-' to read standard input"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero(    "0",    "--rescale must be a positive finite number, but got 0")]
  #[case::negative("-2",   "--rescale must be a positive finite number, but got -2")]
  #[case::infinite("inf",  "--rescale must be a positive finite number, but got inf")]
  #[case::nan(     "NaN",  "--rescale must be a positive finite number, but got NaN")]
  #[trace]
  fn test_args_homoplasy_rejects_invalid_rescale(#[case] value: &str, #[case] expected: &str) -> Result<(), Report> {
    let raw = TreetimeHomoplasyArgsRaw::try_parse_from([
      "homoplasy",
      "--tree=tree.nwk",
      "--alignment=a.fasta",
      &format!("--rescale={value}"),
    ])?;
    assert_error!(TreetimeHomoplasyArgs::try_from(raw), expected);
    Ok(())
  }
}
