#[cfg(test)]
mod tests {
  use crate::partition::marginal::sparse::candidate_sites::SiteSlot;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;

  const SET: bool = true;
  const OR: bool = false;
  const S: SiteSlot = SiteSlot::State;
  const R0: SiteSlot = SiteSlot::Reference(0);
  const R1: SiteSlot = SiteSlot::Reference(1);

  #[rustfmt::skip]
  #[rstest]
  #[case::no_events(                     1, &[],                                                                          &[])]
  #[case::overwrite_after_or_insert(     1, &[(5, S, OR,  b'A'), (5, S, SET, b'C')],                                      &[(5, b'C', "_")])]
  #[case::or_insert_after_overwrite(     1, &[(5, S, SET, b'C'), (5, S, OR,  b'A')],                                      &[(5, b'C', "_")])]
  #[case::first_of_two_or_inserts(       1, &[(5, S, OR,  b'A'), (5, S, OR,  b'G')],                                      &[(5, b'A', "_")])]
  #[case::last_of_two_overwrites(        1, &[(5, S, SET, b'A'), (5, S, SET, b'G')],                                      &[(5, b'G', "_")])]
  #[case::slots_resolve_independently(   1, &[(5, S, OR,  b'A'), (5, R0, SET, b'C'), (5, R0, OR, b'G'), (5, S, SET, b'T')], &[(5, b'T', "C")])]
  #[case::positions_come_out_sorted(     1, &[(9, S, OR,  b'A'), (2, S, OR,  b'C'), (5, S, OR, b'G')],                    &[(2, b'C', "_"), (5, b'G', "_"), (9, b'A', "_")])]
  #[case::interleaved_positions(         1, &[(5, S, OR,  b'A'), (7, S, OR,  b'C'), (5, S, SET, b'G')],                   &[(5, b'G', "_"), (7, b'C', "_")])]
  #[case::reference_only_emits_no_site(  1, &[(3, R0, SET, b'A'), (4, S, OR, b'C')],                                      &[(4, b'C', "_")])]
  #[case::reference_per_message(         2, &[(5, S, SET, b'A'), (5, R1, SET, b'C')],                                     &[(5, b'A', "_C")])]
  #[trace]
  fn test_candidate_sites_events_resolve_as_map_inserts(
    #[case] n_messages: usize,
    #[case] events: &[(usize, SiteSlot, bool, u8)],
    #[case] expected: &[(usize, u8, &str)],
  ) {
    let actual = helpers::summarize(&helpers::finish(n_messages, events));
    assert_eq!(helpers::sites(expected), actual);
  }

  #[test]
  fn test_candidate_sites_fill_missing_references_fills_only_missing() {
    let mut sites = helpers::finish(
      2,
      &[
        (5, S, SET, b'A'),
        (5, R0, SET, b'C'),
        (8, S, OR, b'G'),
        (8, R1, OR, b'T'),
      ],
    );
    sites.fill_missing_references(|message, pos, state| match (message, pos) {
      (0, _) => state,
      (_, 5) => helpers::char(b'-'),
      _ => helpers::char(b'N'),
    });
    assert_eq!(
      helpers::sites(&[(5, b'A', "C-"), (8, b'G', "GT")]),
      helpers::summarize(&sites)
    );
  }

  proptest! {
    #[test]
    fn test_prop_candidate_sites_events_matches_map_model((n_messages, events) in generators::events()) {
      let expected = helpers::map_model(n_messages, &events);
      let actual = helpers::summarize(&helpers::finish(n_messages, &events));
      prop_assert_eq!(expected, actual);
    }
  }

  mod generators {
    use crate::partition::marginal::sparse::candidate_sites::SiteSlot;
    use proptest::collection::vec;
    use proptest::prelude::*;

    pub(super) fn events() -> impl Strategy<Value = (usize, Vec<(usize, SiteSlot, bool, u8)>)> {
      (1..4_usize).prop_flat_map(|n_messages| {
        let event = (
          0..8_usize,
          0..=n_messages,
          any::<bool>(),
          prop::sample::select(b"ACGT".to_vec()),
        )
          .prop_map(|(pos, slot, overwrite, value)| {
            let slot = match slot {
              0 => SiteSlot::State,
              index => SiteSlot::Reference(index - 1),
            };
            (pos, slot, overwrite, value)
          });
        (Just(n_messages), vec(event, 0..40))
      })
    }
  }

  mod helpers {
    use crate::partition::marginal::sparse::candidate_sites::{CandidateSiteEvents, CandidateSites, SiteSlot};
    use std::collections::BTreeMap;
    use treetime_primitives::AsciiChar;

    pub(super) type Site = (usize, AsciiChar, Vec<Option<AsciiChar>>);

    pub(super) fn finish(n_messages: usize, events: &[(usize, SiteSlot, bool, u8)]) -> CandidateSites {
      let mut recorder = CandidateSiteEvents::with_capacity(n_messages, events.len());
      for &(pos, slot, overwrite, value) in events {
        if overwrite {
          recorder.insert(pos, slot, char(value));
        } else {
          recorder.or_insert(pos, slot, char(value));
        }
      }
      recorder.finish()
    }

    pub(super) fn summarize(sites: &CandidateSites) -> Vec<Site> {
      sites
        .iter()
        .map(|(pos, state, references)| (pos, state, references.to_vec()))
        .collect()
    }

    pub(super) fn sites(sites: &[(usize, u8, &str)]) -> Vec<Site> {
      sites
        .iter()
        .map(|&(pos, state, references)| {
          let references = references.bytes().map(|c| (c != b'_').then(|| char(c))).collect();
          (pos, char(state), references)
        })
        .collect()
    }

    pub(super) fn map_model(n_messages: usize, events: &[(usize, SiteSlot, bool, u8)]) -> Vec<Site> {
      let mut slots = vec![BTreeMap::new(); n_messages + 1];
      for &(pos, slot, overwrite, value) in events {
        let index = slot.index();
        let value = char(value);
        if overwrite {
          slots[index].insert(pos, value);
        } else {
          slots[index].entry(pos).or_insert(value);
        }
      }
      slots[0]
        .iter()
        .map(|(&pos, &state)| {
          let references = slots[1..].iter().map(|refs| refs.get(&pos).copied()).collect();
          (pos, state, references)
        })
        .collect()
    }

    pub(super) fn char(c: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(c)
    }
  }
}
