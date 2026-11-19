#[cfg(test)]
mod tests {
  use crate::commands::clock::args::TreetimeClockArgsRaw;
  use clap::Parser;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::default( &[],                                  vec![',', '\t', ';'])]
  #[case::explicit(&["--metadata-delimiters", "|", ":"], vec!['|', ':'])]
  #[trace]
  fn test_metadata_args_delimiters(#[case] extra: &[&str], #[case] expected: Vec<char>) {
    let argv = [&"treetime", &"--metadata=/dev/null"].into_iter().chain(extra);
    let args = TreetimeClockArgsRaw::try_parse_from(argv).unwrap();
    assert_eq!(expected, args.metadata_id.metadata_delimiters);
  }
}
