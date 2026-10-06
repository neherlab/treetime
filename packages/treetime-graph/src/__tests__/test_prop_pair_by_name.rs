#[cfg(test)]
mod tests {
  use crate::node::GraphNodeKey;
  use crate::pair_by_name::pair_by_name;
  use itertools::Itertools;
  use proptest::prelude::*;
  use std::collections::{BTreeMap, BTreeSet};

  proptest! {
    #[test]
    fn test_prop_pair_by_name_pairs_every_named_candidate_with_the_first_entry_of_its_name(
      node_names in prop::collection::vec(prop::option::of("[a-d]"), 0..12),
      candidate_mask in prop::collection::vec(any::<bool>(), 12),
      entry_names in prop::collection::vec("[a-e]", 0..16),
    ) {
      let names: BTreeMap<GraphNodeKey, Option<String>> = node_names
        .iter()
        .enumerate()
        .map(|(index, name)| (GraphNodeKey(index), name.clone()))
        .collect();
      let candidates = names.keys().copied().filter(|key| candidate_mask[key.as_usize()]).collect_vec();
      let entries = entry_names.iter().cloned().enumerate().map(|(index, name)| (name, index)).collect_vec();

      let actual = pair_by_name(candidates.iter().copied(), &names, entries);

      let first_index = |name: &str| entry_names.iter().position(|entry| entry == name);
      let expected_by_node: BTreeMap<GraphNodeKey, usize> = candidates
        .iter()
        .filter_map(|key| Some((*key, first_index(names[key].as_deref()?)?)))
        .collect();
      let candidate_names: BTreeSet<&str> = candidates.iter().filter_map(|key| names[key].as_deref()).collect();
      let expected_unmatched = entry_names
        .iter()
        .unique()
        .filter(|name| !candidate_names.contains(name.as_str()))
        .map(|name| (name.clone(), first_index(name).expect("the name has an entry")))
        .collect_vec();
      let expected_duplicates = entry_names
        .iter()
        .counts()
        .into_iter()
        .filter(|(_, count)| *count > 1)
        .map(|(name, _)| name.clone())
        .sorted()
        .collect_vec();

      prop_assert_eq!(expected_by_node, actual.by_node);
      prop_assert_eq!(expected_unmatched, actual.unmatched);
      prop_assert_eq!(expected_duplicates, actual.duplicate_entry_names);
    }
  }
}
