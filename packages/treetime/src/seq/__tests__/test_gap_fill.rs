#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::TreetimeAncestralArgs;
  use crate::seq::gap_fill::{GapFill, apply_gap_fill};
  use clap::Parser;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_primitives::{AsciiChar, Seq};

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn seq(s: &str) -> Seq {
    Seq::try_from_str(s).unwrap()
  }

  const GAP: u8 = b'-';
  const UNK: u8 = b'N';

  // --- GapFill::OnlyTerminal ---

  #[rustfmt::skip]
  #[rstest]
  #[case::leading_gaps(          "---ACGT",      "NNNACGT")]
  #[case::trailing_gaps(         "ACGT---",      "ACGTNNN")]
  #[case::leading_and_trailing(  "--ACGT--",     "NNACGTNN")]
  #[case::internal_gaps(         "AC--GT",       "AC--GT")]
  #[case::all_gap(               "----",         "NNNN")]
  #[case::no_gap(                "ACGT",         "ACGT")]
  #[case::single_char(           "A",            "A")]
  #[case::single_gap(            "-",            "N")]
  #[case::mixed_terminal_internal("--AC--GT--",  "NNAC--GTNN")]
  #[case::one_non_gap_at_start(  "A---",         "ANNN")]
  #[case::one_non_gap_at_end(    "---A",         "NNNA")]
  #[case::one_non_gap_middle(    "--A--",        "NNANN")]
  #[trace]
  fn test_gap_fill_only_terminal(#[case] input: &str, #[case] expected: &str) {
    let mut s = seq(input);
    apply_gap_fill(&mut s, GapFill::OnlyTerminal, c(GAP), c(UNK));
    assert_eq!(expected, s.as_str());
  }

  #[test]
  fn test_gap_fill_only_terminal_empty() {
    let mut s = seq("");
    apply_gap_fill(&mut s, GapFill::OnlyTerminal, c(GAP), c(UNK));
    assert_eq!("", s.as_str());
  }

  // --- GapFill::All ---

  #[rustfmt::skip]
  #[rstest]
  #[case::leading_gaps(          "---ACGT",      "NNNACGT")]
  #[case::trailing_gaps(         "ACGT---",      "ACGTNNN")]
  #[case::internal_gaps(         "AC--GT",       "ACNNGT")]
  #[case::all_gap(               "----",         "NNNN")]
  #[case::no_gap(                "ACGT",         "ACGT")]
  #[case::single_gap(            "-",            "N")]
  #[case::mixed(                 "--AC--GT--",   "NNACNNGTNN")]
  #[trace]
  fn test_gap_fill_all(#[case] input: &str, #[case] expected: &str) {
    let mut s = seq(input);
    apply_gap_fill(&mut s, GapFill::All, c(GAP), c(UNK));
    assert_eq!(expected, s.as_str());
  }

  #[test]
  fn test_gap_fill_all_empty() {
    let mut s = seq("");
    apply_gap_fill(&mut s, GapFill::All, c(GAP), c(UNK));
    assert_eq!("", s.as_str());
  }

  // --- GapFill::None ---

  #[rustfmt::skip]
  #[rstest]
  #[case::leading_gaps(          "---ACGT",      "---ACGT")]
  #[case::trailing_gaps(         "ACGT---",      "ACGT---")]
  #[case::internal_gaps(         "AC--GT",       "AC--GT")]
  #[case::all_gap(               "----",         "----")]
  #[case::no_gap(                "ACGT",         "ACGT")]
  #[trace]
  fn test_gap_fill_none(#[case] input: &str, #[case] expected: &str) {
    let mut s = seq(input);
    apply_gap_fill(&mut s, GapFill::None, c(GAP), c(UNK));
    assert_eq!(expected, s.as_str());
  }

  // --- Amino acid alphabet (gap='-', unknown='X') ---

  #[rustfmt::skip]
  #[rstest]
  #[case::aa_terminal(       "---ARNDCQ---", "XXXARNDCQXXX")]
  #[case::aa_internal(       "AR--ND",       "AR--ND")]
  #[case::aa_all_gap(        "---",          "XXX")]
  #[trace]
  fn test_gap_fill_amino_acid(#[case] input: &str, #[case] expected: &str) {
    let mut s = seq(input);
    apply_gap_fill(&mut s, GapFill::OnlyTerminal, c(b'-'), c(b'X'));
    assert_eq!(expected, s.as_str());
  }

  // --- v0 parity: OnlyTerminal matches v0 seq2array(fill_overhangs=True) ---

  #[test]
  fn test_gap_fill_v0_parity_typical_sequence() {
    // Simulates incomplete sequencing: leading and trailing gaps with internal variation
    let mut s = seq("----ACGTACGT--ACGT----");
    apply_gap_fill(&mut s, GapFill::OnlyTerminal, c(GAP), c(UNK));
    assert_eq!("NNNNACGTACGT--ACGTNNNN", s.as_str());
  }

  #[test]
  fn test_gap_fill_v0_parity_no_internal_gaps_touched() {
    let mut s = seq("--A-C-G-T--");
    apply_gap_fill(&mut s, GapFill::OnlyTerminal, c(GAP), c(UNK));
    assert_eq!("NNA-C-G-TNN", s.as_str());
  }

  // --- CLI arg resolution ---

  fn base_args() -> Vec<&'static str> {
    vec!["ancestral", "--tree=t.nwk", "--outdir=out"]
  }

  #[test]
  fn test_gap_fill_cli_default_is_only_terminal() {
    let args = TreetimeAncestralArgs::try_parse_from(base_args()).unwrap();
    assert_eq!(GapFill::OnlyTerminal, args.effective_gap_fill());
  }

  #[test]
  fn test_gap_fill_cli_explicit_none() {
    let mut a = base_args();
    a.push("--gap-fill=none");
    let args = TreetimeAncestralArgs::try_parse_from(a).unwrap();
    assert_eq!(GapFill::None, args.effective_gap_fill());
  }

  #[test]
  fn test_gap_fill_cli_explicit_all() {
    let mut a = base_args();
    a.push("--gap-fill=all");
    let args = TreetimeAncestralArgs::try_parse_from(a).unwrap();
    assert_eq!(GapFill::All, args.effective_gap_fill());
  }

  #[test]
  fn test_gap_fill_cli_keep_overhangs_resolves_to_none() {
    let mut a = base_args();
    a.push("--keep-overhangs");
    let args = TreetimeAncestralArgs::try_parse_from(a).unwrap();
    assert_eq!(GapFill::None, args.effective_gap_fill());
  }

  #[test]
  fn test_gap_fill_cli_both_flags_is_error() {
    let mut a = base_args();
    a.push("--keep-overhangs");
    a.push("--gap-fill=all");
    let result = TreetimeAncestralArgs::try_parse_from(a);
    result.unwrap_err();
  }
}
