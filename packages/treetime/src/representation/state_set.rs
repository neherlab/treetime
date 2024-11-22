use serde::{Deserialize, Serialize};
use std::borrow::Borrow;
use std::collections::btree_set::Iter;
use std::collections::BTreeSet;

#[derive(Clone, Debug)]
pub enum StateSetStatus<'a> {
  Empty,
  Unambiguous(char),
  Ambiguous(Iter<'a, char>),
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct StateSet {
  data: BTreeSet<char>,
}

impl StateSet {
  /// Create an empty set
  pub fn new() -> Self {
    Self { data: BTreeSet::new() }
  }

  /// Create a set containing a given chars
  pub fn from_chars<I, T>(chars: I) -> Self
  where
    I: IntoIterator<Item = T>,
    T: Borrow<char>,
  {
    Self {
      data: chars.into_iter().map(|c| *c.borrow()).collect(),
    }
  }

  /// Create a set containing a single given char
  pub fn from_char(c: char) -> Self {
    Self::from_chars([c])
  }

  /// Create a set from intersection of given sets
  pub fn from_intersection(sets: &[StateSet]) -> StateSet {
    let intersection = sets
      .iter()
      .map(|set| &set.data)
      .fold(None, |acc: Option<BTreeSet<char>>, set| {
        if let Some(a) = acc {
          Some(a.intersection(set).copied().collect())
        } else {
          Some(set.clone())
        }
      })
      .unwrap_or_else(BTreeSet::new);

    StateSet { data: intersection }
  }

  /// Create a set from union of given sets
  pub fn from_union(sets: &[StateSet]) -> StateSet {
    let data = sets
      .iter()
      .flat_map(|set| set.data.iter())
      .copied()
      .collect::<BTreeSet<char>>();
    StateSet { data }
  }

  /// Get contents - empty, unique element or multiple ambiguous elements
  pub fn get(&self) -> StateSetStatus {
    match self.data.len() {
      0 => StateSetStatus::Empty,
      1 => StateSetStatus::Unambiguous(*self.data.iter().next().unwrap()),
      _ => StateSetStatus::Ambiguous(self.data.iter()),
    }
  }

  /// Get one of the characters from the set, or panic if the set is empty.
  /// Note: which of the multiple character gets retrieved is not specified and is not to be relied on.
  pub fn get_one(&self) -> char {
    assert!(self.data.len() > 0);
    self.data.iter().next().copied().unwrap()
  }

  /// Get exactly one element from the set, or panic if there's not exactly one element.
  /// Note: which of the multiple character gets retrieved is not specified and is not to be relied on.
  pub fn get_one_exactly(&self) -> char {
    assert_eq!(self.data.len(), 1);
    self.data.iter().next().copied().unwrap()
  }

  /// Check if set contains a character
  pub fn contains(&self, c: char) -> bool {
    self.data.contains(&c)
  }

  /// Access characters
  #[allow(clippy::needless_lifetimes)]
  pub fn chars<'a>(&'a self) -> impl Iterator<Item = char> + 'a {
    self.data.iter().copied()
  }

  /// Access internal set implementation
  pub fn inner(&self) -> BTreeSet<char> {
    self.data.clone()
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use maplit::btreeset;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_stateset_new() {
    let set = StateSet::new();
    let actual = set.inner();
    let expected: BTreeSet<char> = BTreeSet::new();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_stateset_from_chars() {
    let chars = btreeset!['A', 'C', 'G', 'T'];
    let set = StateSet::from_chars(&chars);
    let actual = &set.inner();
    let expected = &chars;
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_stateset_from_char() {
    let set = StateSet::from_char('A');
    let actual = set.inner();
    let expected = btreeset! {'A'};
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_stateset_from_intersection() {
    let set1 = StateSet::from_chars(vec!['A', 'C', 'G']);
    let set2 = StateSet::from_chars(vec!['C', 'G', 'T']);
    let result = StateSet::from_intersection(&[set1, set2]);
    let actual = result.inner();
    let expected = btreeset! {'C', 'G'};
    assert_eq!(expected, actual);

    let empty_intersection = StateSet::from_intersection(&[]);
    let actual = empty_intersection.inner();
    let expected: BTreeSet<char> = BTreeSet::new();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_stateset_from_union() {
    let set1 = StateSet::from_chars(vec!['A', 'C']);
    let set2 = StateSet::from_chars(vec!['G', 'T']);
    let result = StateSet::from_union(&[set1, set2]);
    let actual = result.inner();
    let expected = btreeset! {'A', 'C', 'G', 'T'};
    assert_eq!(expected, actual);

    let empty_union = StateSet::from_union(&[]);
    let actual = empty_union.inner();
    let expected: BTreeSet<char> = BTreeSet::new();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_stateset_get_empty() {
    let set = StateSet::new();
    assert!(matches!(set.get(), StateSetStatus::Empty));
  }

  #[test]
  fn test_stateset_get_unambiguous() {
    let set = StateSet::from_char('A');
    assert!(matches!(set.get(), StateSetStatus::Unambiguous('A')));
  }

  #[test]
  fn test_stateset_get_ambiguous() {
    let set = StateSet::from_chars(vec!['A', 'C']);
    assert!(matches!(set.get(), StateSetStatus::Ambiguous(_)));
    if let StateSetStatus::Ambiguous(iter) = set.get() {
      let actual: BTreeSet<char> = iter.copied().collect();
      let expected = btreeset! {'A', 'C'};
      assert_eq!(expected, actual);
    }
  }

  #[test]
  fn test_stateset_get_one() {
    let set = StateSet::from_chars(['T', 'A']);
    let actual = set.get_one();
    let expected = 'A';
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_stateset_get_one_exactly() {
    let set = StateSet::from_chars(['T']);
    let actual = set.get_one();
    let expected = 'T';
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_stateset_inner() {
    let chars = vec!['A', 'G', 'T'];
    let set = StateSet::from_chars(chars.clone());
    let actual = set.inner();
    let expected = chars.into_iter().collect::<BTreeSet<char>>();
    assert_eq!(expected, actual);
  }
}
