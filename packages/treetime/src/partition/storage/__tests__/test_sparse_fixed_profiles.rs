#[cfg(test)]
mod tests {
  use crate::partition::storage::sparse::FixedProfiles;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use smallvec::smallvec;
  use std::collections::BTreeMap;

  #[test]
  fn test_sparse_fixed_profiles_iterate_in_state_order_and_replace_on_insert() {
    let mut fixed = FixedProfiles::default();
    fixed.insert(helpers::char(b'G'), 3, smallvec![0.0, 0.0, 1.0, 0.0]);
    fixed.insert(helpers::char(b'A'), 5, smallvec![1.0, 0.0, 0.0, 0.0]);
    fixed.insert(helpers::char(b'T'), 1, smallvec![0.0, 0.0, 0.0, 1.0]);
    fixed.insert(helpers::char(b'A'), 2, smallvec![0.5, 0.5, 0.0, 0.0]);

    let expected = vec![
      (helpers::char(b'A'), 2, vec![0.5, 0.5, 0.0, 0.0]),
      (helpers::char(b'G'), 3, vec![0.0, 0.0, 1.0, 0.0]),
      (helpers::char(b'T'), 1, vec![0.0, 0.0, 0.0, 1.0]),
    ];
    assert_eq!(expected, helpers::entries(&fixed));
  }

  proptest! {
    #[test]
    fn test_prop_sparse_fixed_profiles_matches_btreemap(inserts in generators::inserts()) {
      let mut fixed = FixedProfiles::default();
      let mut model = BTreeMap::new();
      for (state, count, dis) in inserts {
        fixed.insert(state, count, dis.iter().copied().collect());
        model.insert(state, (count, dis));
      }

      let expected_lookups = generators::STATES.iter().map(|&c| model.get(&helpers::char(c)).cloned()).collect_vec();
      let actual_lookups = generators::STATES
        .iter()
        .map(|&c| fixed.get(helpers::char(c)).map(|entry| (entry.count, entry.dis.to_vec())))
        .collect_vec();
      prop_assert_eq!(expected_lookups, actual_lookups);

      let expected = model.into_iter().map(|(state, (count, dis))| (state, count, dis)).collect_vec();
      prop_assert_eq!(expected, helpers::entries(&fixed));
    }
  }

  mod generators {
    use proptest::collection::vec;
    use proptest::prelude::*;
    use treetime_primitives::AsciiChar;

    pub(super) const STATES: &[u8] = b"-ACGNRTY";

    pub(super) fn inserts() -> impl Strategy<Value = Vec<(AsciiChar, usize, Vec<f64>)>> {
      let insert = (
        prop::sample::select(STATES.to_vec()).prop_map(AsciiChar::from_byte_unchecked),
        0..100_usize,
        vec(0.0..1.0_f64, 4),
      );
      vec(insert, 0..20)
    }
  }

  mod helpers {
    use crate::partition::storage::sparse::FixedProfiles;
    use treetime_primitives::AsciiChar;

    pub(super) fn entries(fixed: &FixedProfiles) -> Vec<(AsciiChar, usize, Vec<f64>)> {
      fixed
        .iter()
        .map(|entry| (entry.state, entry.count, entry.dis.to_vec()))
        .collect()
    }

    pub(super) fn char(c: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(c)
    }
  }
}
