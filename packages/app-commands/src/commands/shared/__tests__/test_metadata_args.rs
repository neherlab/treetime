#[cfg(test)]
mod tests {
  use crate::commands::clock::args::TreetimeClockArgsRaw;
  use clap::Parser;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_metadata_args_default_delimiters() {
    let args = TreetimeClockArgsRaw::try_parse_from(["treetime", "--metadata=/dev/null"]).unwrap();
    assert_eq!(vec![',', '\t', ';'], args.metadata_id.metadata_delimiters);
  }

  #[test]
  fn test_metadata_args_explicit_delimiters() {
    let args =
      TreetimeClockArgsRaw::try_parse_from(["treetime", "--metadata=/dev/null", "--metadata-delimiters", "|", ":"])
        .unwrap();
    assert_eq!(vec!['|', ':'], args.metadata_id.metadata_delimiters);
  }
}
