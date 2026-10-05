#[cfg(test)]
mod tests {
  use crate::commands::homoplasy::args::{TreetimeHomoplasyArgs, TreetimeHomoplasyArgsRaw};
  use clap::Parser;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[test]
  fn test_args_homoplasy_reads_homoplasy_flags() -> Result<(), Report> {
    let raw = TreetimeHomoplasyArgsRaw::try_parse_from([
      "homoplasy",
      "--tree=tree.nwk",
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
        args.zero_based
      )
    );
    Ok(())
  }

  #[test]
  fn test_args_homoplasy_defaults() -> Result<(), Report> {
    let raw = TreetimeHomoplasyArgsRaw::try_parse_from(["homoplasy", "--tree=tree.nwk"])?;
    let args = TreetimeHomoplasyArgs::try_from(raw)?;

    assert_eq!(
      (0, 1.0, false, None, 10, false),
      (
        args.constant_sites,
        args.rescale,
        args.detailed,
        args.drms,
        args.num_mut,
        args.zero_based
      )
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
    let raw = TreetimeHomoplasyArgsRaw::try_parse_from(["homoplasy", "--tree=tree.nwk", &format!("--rescale={value}")])?;
    assert_error!(TreetimeHomoplasyArgs::try_from(raw), expected);
    Ok(())
  }
}
