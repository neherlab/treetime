#[cfg(test)]
mod tests {
  use crate::seq::overlay::SeqOverlay;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;
  use std::sync::Arc;
  use treetime_primitives::Seq;

  #[test]
  fn test_overlay_to_seq_replaces_states_at_overridden_positions() -> Result<(), Report> {
    let overlay = helpers::overlay("ACGTAC", "T...G.")?;
    assert_eq!(Seq::try_from_str("TCGTGC")?, overlay.to_seq());
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::equal_without_overrides(  ("ACGT",   "...."),   ("ACGT",   "...."),   0)]
  #[case::same_override(            ("ACGT",   ".T.."),   ("ACGT",   ".T.."),   0)]
  #[case::different_override(       ("ACGT",   ".T.."),   ("ACGT",   ".G.."),   1)]
  #[case::override_equal_to_base(   ("ACGT",   ".C.."),   ("ACGT",   "...."),   0)]
  #[case::disjoint_overrides(       ("ACGT",   "G..."),   ("ACGT",   "...A"),   2)]
  #[case::different_bases(          ("ACGT",   "...."),   ("TCGA",   "...."),   2)]
  #[case::override_hides_base_diff( ("ACGT",   "T..."),   ("TCGT",   "...."),   0)]
  #[case::different_lengths(        ("AAAA",   "...."),   ("AAAAAA", "......"), 2)]
  #[trace]
  fn test_overlay_count_differences_counts_differing_and_extra_positions(
    #[case] (base_a, mask_a): (&str, &str),
    #[case] (base_b, mask_b): (&str, &str),
    #[case] expected: usize,
  ) -> Result<(), Report> {
    let a = helpers::overlay(base_a, mask_a)?;
    let b = helpers::overlay(base_b, mask_b)?;
    assert_eq!(expected, a.count_differences(&b));
    assert_eq!(expected, b.count_differences(&a));
    Ok(())
  }

  #[test]
  fn test_overlay_count_differences_of_a_shared_base_compares_overrides() -> Result<(), Report> {
    let base = Arc::new(Seq::try_from_str("ACGTACGT")?);
    let a = SeqOverlay::new(Arc::clone(&base), helpers::overrides("..A..T.."));
    let b = SeqOverlay::new(base, helpers::overrides("..A....C"));
    assert_eq!(2, a.count_differences(&b));
    Ok(())
  }

  proptest! {
    #[test]
    fn test_prop_overlay_count_differences_equals_count_of_materialized_sequences(
      (base, overrides_a, overrides_b, other_base) in generators::overlay_pair(),
    ) {
      let base = Arc::new(base);
      let a = SeqOverlay::new(Arc::clone(&base), overrides_a);
      let b = SeqOverlay::new(other_base.map_or(base, Arc::new), overrides_b);
      let (seq_a, seq_b) = (a.to_seq(), b.to_seq());
      let expected = seq_a.iter().zip(seq_b.iter()).filter(|(x, y)| x != y).count() + seq_a.len().abs_diff(seq_b.len());
      prop_assert_eq!(expected, a.count_differences(&b));
    }
  }

  mod generators {
    use proptest::collection::{btree_map, vec};
    use proptest::prelude::*;
    use std::collections::BTreeMap;
    use treetime_primitives::{AsciiChar, Seq};

    pub(super) fn overlay_pair()
    -> impl Strategy<Value = (Seq, BTreeMap<usize, AsciiChar>, BTreeMap<usize, AsciiChar>, Option<Seq>)> {
      (1_usize..40).prop_flat_map(|length| {
        (
          seq(length),
          overrides(length),
          overrides(length),
          prop::option::of(seq(length)),
        )
      })
    }

    fn seq(length: usize) -> impl Strategy<Value = Seq> {
      vec(state(), length).prop_map(|states| states.into_iter().collect())
    }

    fn overrides(length: usize) -> impl Strategy<Value = BTreeMap<usize, AsciiChar>> {
      btree_map(0..length, state(), 0..length)
    }

    fn state() -> impl Strategy<Value = AsciiChar> {
      prop::sample::select(b"ACGT".to_vec()).prop_map(AsciiChar::from_byte_unchecked)
    }
  }

  mod helpers {
    use crate::seq::overlay::SeqOverlay;
    use eyre::Report;
    use std::collections::BTreeMap;
    use std::sync::Arc;
    use treetime_primitives::{AsciiChar, Seq};

    pub(super) fn overlay(base: &str, mask: &str) -> Result<SeqOverlay, Report> {
      Ok(SeqOverlay::new(Arc::new(Seq::try_from_str(base)?), overrides(mask)))
    }

    pub(super) fn overrides(mask: &str) -> BTreeMap<usize, AsciiChar> {
      mask
        .bytes()
        .enumerate()
        .filter(|&(_, state)| state != b'.')
        .map(|(pos, state)| (pos, AsciiChar::from_byte_unchecked(state)))
        .collect()
    }
  }
}
