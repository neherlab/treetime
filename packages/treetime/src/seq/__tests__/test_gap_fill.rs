#[cfg(test)]
mod tests {
  use crate::seq::gap_fill::{GapFill, apply_gap_fill};
  use helpers::{c, seq};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::terminal_empty(               GapFill::OnlyTerminal, "",           "")]
  #[case::terminal_no_gap(              GapFill::OnlyTerminal, "ACGT",       "ACGT")]
  #[case::terminal_leading(             GapFill::OnlyTerminal, "---ACGT",    "NNNACGT")]
  #[case::terminal_trailing(            GapFill::OnlyTerminal, "ACGT---",    "ACGTNNN")]
  #[case::terminal_keeps_internal(      GapFill::OnlyTerminal, "AC--GT",     "AC--GT")]
  #[case::terminal_all_gap(             GapFill::OnlyTerminal, "----",       "NNNN")]
  #[case::terminal_single_gap(          GapFill::OnlyTerminal, "-",          "N")]
  #[case::terminal_mixed(               GapFill::OnlyTerminal, "--AC--GT--", "NNAC--GTNN")]
  #[case::terminal_alternating(         GapFill::OnlyTerminal, "--A-C-G-T--", "NNA-C-G-TNN")]
  #[case::terminal_one_non_gap_at_start(GapFill::OnlyTerminal, "A---",       "ANNN")]
  #[case::terminal_one_non_gap_at_end(  GapFill::OnlyTerminal, "---A",       "NNNA")]
  #[case::terminal_one_non_gap_middle(  GapFill::OnlyTerminal, "--A--",      "NNANN")]
  #[case::all_empty(                    GapFill::All,          "",           "")]
  #[case::all_no_gap(                   GapFill::All,          "ACGT",       "ACGT")]
  #[case::all_mixed(                    GapFill::All,          "--AC--GT--", "NNACNNGTNN")]
  #[case::all_all_gap(                  GapFill::All,          "----",       "NNNN")]
  #[case::none_keeps_every_gap(         GapFill::None,         "--AC--GT--", "--AC--GT--")]
  #[trace]
  fn test_gap_fill_replaces_the_gaps_that_the_mode_selects(
    #[case] mode: GapFill,
    #[case] input: &str,
    #[case] expected: &str,
  ) {
    let mut actual = seq(input);
    apply_gap_fill(&mut actual, mode, c(b'-'), c(b'N'));
    assert_eq!(expected, actual.as_str());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::terminal_gaps(  "---ARNDCQ---", "XXXARNDCQXXX")]
  #[case::all_gap(        "---",          "XXX")]
  #[trace]
  fn test_gap_fill_writes_the_given_unknown_character(#[case] input: &str, #[case] expected: &str) {
    let mut actual = seq(input);
    apply_gap_fill(&mut actual, GapFill::OnlyTerminal, c(b'-'), c(b'X'));
    assert_eq!(expected, actual.as_str());
  }

  mod helpers {
    use treetime_primitives::{AsciiChar, Seq};

    pub(super) fn c(b: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(b)
    }

    pub(super) fn seq(s: &str) -> Seq {
      Seq::try_from_str(s).unwrap()
    }
  }
}
