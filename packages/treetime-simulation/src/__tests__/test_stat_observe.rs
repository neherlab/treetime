#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::helpers::{chisquare_p_value, three_tip_tree};
  use crate::FALSE_FAIL_RATE;
  use crate::observe::{TipObservation, codes_containing, observe_tips};
  use itertools::Itertools;
  use rstest::rstest;
  use std::iter::repeat_n;
  use treetime_primitives::{AsciiChar, Seq};
  use treetime_utils::sync::random::get_random_number_generator;

  const SITE_COUNT: usize = 60_000;

  #[rustfmt::skip]
  #[rstest]
  #[case::a(b'A')]
  #[case::c(b'C')]
  #[case::g(b'G')]
  #[case::t(b'T')]
  #[trace]
  fn test_stat_observe_tips_draw_each_code_with_the_true_base_equally_often(#[case] base: u8) {
    let tree = three_tip_tree();
    let base = AsciiChar::from_byte_unchecked(base);
    let seq = repeat_n(base, SITE_COUNT).collect::<Seq>();
    let sequences = [tree.tip_a, tree.tip_b, tree.tip_c].into_iter().map(|leaf| (leaf, seq.clone())).collect();
    let observation = TipObservation {
      ambiguity_rate: 1.0,
      unknown_run_rate: 0.0,
      unknown_run_lengths: 1..=1,
    };

    let observed = observe_tips(&tree.graph, &sequences, &observation, &mut get_random_number_generator(50)).unwrap();

    let codes = codes_containing(base).unwrap();
    let counts = observed.values().flat_map(|seq| seq.iter().copied()).counts();
    let code_counts = codes.iter().map(|code| counts.get(code).copied().unwrap_or(0)).collect_vec();
    let p_value = chisquare_p_value(&code_counts, &vec![1.0 / 6.0; codes.len()]);
    assert!(
      p_value > FALSE_FAIL_RATE,
      "chi-square test of the codes of {SITE_COUNT} sites of {base} against equal probabilities: counts {code_counts:?}, p = {p_value}"
    );
  }
}
