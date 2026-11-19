#[cfg(test)]
mod tests {
  use crate::commands::shared::gap_fill::GapFillArgs;
  use clap::Parser;
  use helpers::effective;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::seq::gap_fill::GapFill;

  #[derive(Debug, Parser)]
  struct GapFillCli {
    #[command(flatten)]
    gap_fill_args: GapFillArgs,
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::default_is_only_terminal(&[],                  GapFill::OnlyTerminal)]
  #[case::explicit_all(            &["--gap-fill=all"],  GapFill::All)]
  #[case::explicit_none(           &["--gap-fill=none"], GapFill::None)]
  #[trace]
  fn test_gap_fill_args_command_line_selects_the_policy(#[case] extra: &[&str], #[case] expected: GapFill) {
    assert_eq!(expected, effective(extra));
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
