#[cfg(test)]
mod tests {
  use crate::commands::homoplasy::args::{TreetimeHomoplasyArgs, TreetimeHomoplasyArgsRaw};
  use clap::Parser;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::seq::mutation::PositionBase;
  use treetime_utils::assert_error;

  #[test]
  fn test_args_homoplasy_reads_homoplasy_flags() -> Result<(), Report> {
    let raw = TreetimeHomoplasyArgsRaw::try_parse_from([
      "homoplasy",
      "--tree=tree.nwk",
      "--alignment=a.fasta",
      "--const=100",
      "--rescale=0.5",
      "--detailed",
      "--drms=drms.tsv",
      "-n",
      "3",
      "--zero-based",
    ])?;
    let args = TreetimeHomoplasyArgs::try_from(raw)?;

    assert_eq!(
      (100, 0.5, true, Some("drms.tsv"), 3, true),
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
  fn test_args_homoplasy_defaults() -> Result<(), Report> {
    let raw = TreetimeHomoplasyArgsRaw::try_parse_from(["homoplasy", "--tree=tree.nwk", "--alignment=a.fasta"])?;
    let args = TreetimeHomoplasyArgs::try_from(raw)?;

    assert_eq!(
      (0, 1.0, false, None, 10, false),
      (
        args.constant_sites,
        args.rescale,
        args.detailed,
        args.drms.clone(),
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
