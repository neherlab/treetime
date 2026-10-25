#[cfg(test)]
mod tests {
  use self::helpers::{as_btree_map, key};
  use crate::key_map::NodeMap;
  use crate::keymap;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_key_map_iterates_in_key_order_from_both_ends() {
    let map: NodeMap<&str> = keymap! { key(4) => "d", key(0) => "a", key(2) => "c" };

    let forward = map.iter().map(|(key, &value)| (key, value)).collect_vec();
    let backward = map.iter().rev().map(|(key, &value)| (key, value)).collect_vec();

    assert_eq!(vec![(key(0), "a"), (key(2), "c"), (key(4), "d")], forward);
    assert_eq!(vec![(key(4), "d"), (key(2), "c"), (key(0), "a")], backward);
  }

  #[test]
  fn test_key_map_equality_ignores_removed_entries() {
    let mut with_removed: NodeMap<u8> = keymap! { key(1) => 10, key(7) => 70 };
    with_removed.remove(key(7));
    let without: NodeMap<u8> = keymap! { key(1) => 10 };

    assert_eq!(without, with_removed);
  }

  #[test]
  fn test_key_map_debug_lists_entries_as_a_map() {
    let mut map: NodeMap<u8> = keymap! { key(1) => 10, key(3) => 30 };
    map.remove(key(3));

    assert_eq!("{GraphNodeKey(1): 10}", format!("{map:?}"));
  }

  #[test]
  fn test_key_map_entry_keeps_an_existing_value() {
    let mut map: NodeMap<u8> = keymap! { key(2) => 20 };

    *map.entry(key(2)).or_insert(99) += 1;
    *map.entry(key(5)).or_insert(50) += 1;

    assert_eq!(btreemap! { key(2) => 21, key(5) => 51 }, as_btree_map(&map));
  }

  #[test]
  #[should_panic(expected = "Node 3 has no entry in the map")]
  fn test_key_map_index_panics_for_a_key_without_entry() {
    let map: NodeMap<u8> = keymap! { key(1) => 10 };

    let _value = map[key(3)];
  }

  mod helpers {
    use crate::key_map::NodeMap;
    use crate::node::GraphNodeKey;
    use std::collections::BTreeMap;

    pub(super) const fn key(index: usize) -> GraphNodeKey {
      GraphNodeKey(index)
    }

    pub(super) fn as_btree_map<V: Clone>(map: &NodeMap<V>) -> BTreeMap<GraphNodeKey, V> {
      map.iter().map(|(key, value)| (key, value.clone())).collect()
    }
  }
}
