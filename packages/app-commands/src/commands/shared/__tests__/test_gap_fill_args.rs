#[cfg(test)]
mod tests {
  use crate::commands::shared::gap_fill::GapFillArgs;
  use clap::Parser;
  use helpers::effective;
  use pretty_assertions::assert_eq;
  use treetime::seq::gap_fill::GapFill;

  #[derive(Debug, Parser)]
  struct GapFillCli {
    #[command(flatten)]
    gap_fill_args: GapFillArgs,
  }

  #[test]
  fn test_gap_fill_args_default_is_only_terminal() {
    assert_eq!(GapFill::OnlyTerminal, effective(&[]));
  }

  #[test]
  fn test_gap_fill_args_explicit_none() {
    assert_eq!(GapFill::None, effective(&["--gap-fill=none"]));
  }

  #[test]
  fn test_gap_fill_args_explicit_all() {
    assert_eq!(GapFill::All, effective(&["--gap-fill=all"]));
  }

  mod helpers {
    use super::*;

    pub(super) fn effective(extra: &[&str]) -> GapFill {
      let mut argv = vec!["treetime"];
      argv.extend_from_slice(extra);
      GapFill::from(GapFillCli::try_parse_from(argv).unwrap().gap_fill_args.gap_fill)
    }
  }
}
