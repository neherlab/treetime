#[cfg(test)]
mod tests {
  use crate::partition::storage::var_pos_map::VarPosMap;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use std::collections::BTreeMap;

  #[test]
  fn test_var_pos_map_iterates_by_position_and_replaces_on_insert() {
    let mut map = VarPosMap::default();
    map.insert(7, helpers::char(b'G'), &[0.0, 0.0, 1.0, 0.0]);
    map.insert(2, helpers::char(b'A'), &[1.0, 0.0, 0.0, 0.0]);
    map.insert_with(9, helpers::char(b'T'), 4, |row| row[3] = 1.0);
    map.insert(2, helpers::char(b'C'), &[0.0, 1.0, 0.0, 0.0]);

    let expected = vec![
      (2, helpers::char(b'C'), vec![0.0, 1.0, 0.0, 0.0]),
      (7, helpers::char(b'G'), vec![0.0, 0.0, 1.0, 0.0]),
      (9, helpers::char(b'T'), vec![0.0, 0.0, 0.0, 1.0]),
    ];
    assert_eq!(expected, helpers::entries(&map));
  }

  proptest! {
    #[test]
    fn test_prop_var_pos_map_matches_btreemap(inserts in generators::inserts()) {
      let mut map = VarPosMap::default();
      let mut model = BTreeMap::new();
      for (pos, state, dis) in inserts {
        map.insert(pos, state, &dis);
        model.insert(pos, (state, dis));
      }

      let expected_lookups = (0..generators::N_POSITIONS).map(|pos| model.get(&pos).cloned()).collect_vec();
      let actual_lookups = (0..generators::N_POSITIONS)
        .map(|pos| map.get(pos).map(|var| (var.state, var.dis.to_vec())))
        .collect_vec();
      prop_assert_eq!(expected_lookups, actual_lookups);
      let expected_keys = (0..generators::N_POSITIONS).map(|pos| model.contains_key(&pos)).collect_vec();
      let actual_keys = (0..generators::N_POSITIONS).map(|pos| map.contains_key(pos)).collect_vec();
      prop_assert_eq!(expected_keys, actual_keys);
      prop_assert_eq!(model.len(), map.len());

      let expected = model.into_iter().map(|(pos, (state, dis))| (pos, state, dis)).collect_vec();
      prop_assert_eq!(expected, helpers::entries(&map));
    }
  }

  mod generators {
    use proptest::collection::vec;
    use proptest::prelude::*;
    use treetime_primitives::AsciiChar;

    pub(super) const N_POSITIONS: usize = 20;

    pub(super) fn inserts() -> impl Strategy<Value = Vec<(usize, AsciiChar, Vec<f64>)>> {
      let insert = (
        0..N_POSITIONS,
        prop::sample::select(b"ACGT".to_vec()).prop_map(AsciiChar::from_byte_unchecked),
        vec(0.0..1.0_f64, 4),
      );
      vec(insert, 0..40)
    }
  }

  mod helpers {
    use crate::partition::storage::var_pos_map::VarPosMap;
    use treetime_primitives::AsciiChar;

    pub(super) fn entries(map: &VarPosMap) -> Vec<(usize, AsciiChar, Vec<f64>)> {
      map
        .iter()
        .map(|(&pos, var)| (pos, var.state, var.dis.to_vec()))
        .collect()
    }

    pub(super) fn char(c: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(c)
    }
  }
}
