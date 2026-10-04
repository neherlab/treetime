use crate::partition::storage::sparse::VarPos;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::iter::Map;
use std::ops::Index;
use std::slice::Iter;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(from = "BTreeMap<usize, VarPos>", into = "BTreeMap<usize, VarPos>")]
pub struct VarPosMap {
  entries: Vec<(usize, VarPos)>,
}

impl VarPosMap {
  pub(crate) fn get(&self, pos: usize) -> Option<&VarPos> {
    self
      .entries
      .binary_search_by_key(&pos, |(key, _)| *key)
      .ok()
      .map(|index| &self.entries[index].1)
  }

  pub(crate) fn contains_key(&self, pos: usize) -> bool {
    self.entries.binary_search_by_key(&pos, |(key, _)| *key).is_ok()
  }

  pub(crate) fn insert(&mut self, pos: usize, value: VarPos) {
    if self.entries.last().is_none_or(|(last, _)| *last < pos) {
      self.entries.push((pos, value));
      return;
    }
    match self.entries.binary_search_by_key(&pos, |(key, _)| *key) {
      Ok(index) => self.entries[index].1 = value,
      Err(index) => self.entries.insert(index, (pos, value)),
    }
  }

  pub(crate) fn iter(&self) -> impl Iterator<Item = (&usize, &VarPos)> {
    self.entries.iter().map(|(pos, value)| (pos, value))
  }

  pub(crate) fn keys(&self) -> impl Iterator<Item = &usize> {
    self.entries.iter().map(|(pos, _)| pos)
  }

  pub(crate) fn values(&self) -> impl Iterator<Item = &VarPos> {
    self.entries.iter().map(|(_, value)| value)
  }
}

impl FromIterator<(usize, VarPos)> for VarPosMap {
  fn from_iter<I: IntoIterator<Item = (usize, VarPos)>>(iter: I) -> Self {
    let mut map = Self::default();
    for (pos, value) in iter {
      map.insert(pos, value);
    }
    map
  }
}

impl<'a> IntoIterator for &'a VarPosMap {
  type Item = (&'a usize, &'a VarPos);
  type IntoIter = Map<Iter<'a, (usize, VarPos)>, fn(&'a (usize, VarPos)) -> (&'a usize, &'a VarPos)>;

  fn into_iter(self) -> Self::IntoIter {
    self.entries.iter().map(|(pos, value)| (pos, value))
  }
}

impl From<BTreeMap<usize, VarPos>> for VarPosMap {
  fn from(map: BTreeMap<usize, VarPos>) -> Self {
    Self {
      entries: map.into_iter().collect(),
    }
  }
}

impl From<VarPosMap> for BTreeMap<usize, VarPos> {
  fn from(map: VarPosMap) -> Self {
    map.entries.into_iter().collect()
  }
}

impl Index<usize> for VarPosMap {
  type Output = VarPos;

  #[allow(
    clippy::expect_used,
    reason = "indexing by a missing position is a caller bug, as with BTreeMap"
  )]
  fn index(&self, pos: usize) -> &VarPos {
    self.get(pos).expect("position is not a variable position")
  }
}
