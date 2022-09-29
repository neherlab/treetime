use itertools::Itertools;
use std::collections::HashMap;
use std::hash::Hash;

/// Count occurrences (multiplicities) of unique values in an iterator
pub fn count_occurences<T: Copy + Hash + Eq + Ord>(it: impl Iterator<Item = T>) -> Vec<(T, usize)> {
  let mut occurences: HashMap<T, usize> = HashMap::new();
  for x in it {
    *occurences.entry(x).or_default() += 1;
  }
  occurences.into_iter().sorted_by_key(|(key, _)| *key).collect_vec()
}
