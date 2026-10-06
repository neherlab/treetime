use crate::node::GraphNodeKey;
use indexmap::IndexMap;
use indexmap::map::Entry;
use itertools::izip;
use std::collections::{BTreeMap, BTreeSet};

pub fn pair_by_name<T: Clone>(
  candidates: impl IntoIterator<Item = GraphNodeKey>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  entries: impl IntoIterator<Item = (String, T)>,
) -> NamePairing<T> {
  let mut duplicate_entry_names = BTreeSet::new();
  let mut by_name: IndexMap<String, T> = IndexMap::new();
  for (name, value) in entries {
    match by_name.entry(name) {
      Entry::Occupied(entry) => {
        duplicate_entry_names.insert(entry.key().clone());
      },
      Entry::Vacant(entry) => {
        entry.insert(value);
      },
    }
  }

  let mut keys_by_entry: Vec<Vec<GraphNodeKey>> = vec![vec![]; by_name.len()];
  for key in candidates {
    if let Some(index) = names[&key].as_deref().and_then(|name| by_name.get_index_of(name)) {
      keys_by_entry[index].push(key);
    }
  }

  let mut by_node = BTreeMap::new();
  let mut unmatched = vec![];
  for ((name, value), keys) in izip!(by_name, keys_by_entry) {
    let Some((last, others)) = keys.split_last() else {
      unmatched.push((name, value));
      continue;
    };
    by_node.extend(others.iter().map(|key| (*key, value.clone())));
    by_node.insert(*last, value);
  }

  NamePairing {
    by_node,
    unmatched,
    duplicate_entry_names: duplicate_entry_names.into_iter().collect(),
  }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NamePairing<T> {
  pub by_node: BTreeMap<GraphNodeKey, T>,
  pub unmatched: Vec<(String, T)>,
  pub duplicate_entry_names: Vec<String>,
}
