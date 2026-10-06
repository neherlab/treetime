use crate::node::GraphNodeKey;
use indexmap::IndexMap;
use indexmap::map::Entry;
use itertools::Itertools;
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

  let mut matches = candidates
    .into_iter()
    .filter_map(|key| Some((by_name.get_index_of(names[&key].as_deref()?)?, key)))
    .collect_vec();
  matches.sort_by_key(|(index, _)| *index);
  let mut matches = matches.into_iter().peekable();

  let mut by_node = BTreeMap::new();
  let mut unmatched = vec![];
  for (index, (name, value)) in by_name.into_iter().enumerate() {
    let keys = matches
      .peeking_take_while(|(match_index, _)| *match_index == index)
      .map(|(_, key)| key)
      .collect_vec();
    let Some((last, others)) = keys.split_last() else {
      unmatched.push(name);
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
  pub unmatched: Vec<String>,
  pub duplicate_entry_names: Vec<String>,
}
