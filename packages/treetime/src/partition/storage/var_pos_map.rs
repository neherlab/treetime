use crate::partition::storage::sparse::{SiteProfile, VarPos};
use smallvec::smallvec;
use std::iter::Map;
use std::slice::Iter;
use treetime_primitives::AsciiChar;

impl<'a> IntoIterator for &'a VarPosMap {
  type Item = (&'a usize, VarPos<'a>);
  type IntoIter = Map<Iter<'a, Entry>, fn(&'a Entry) -> (&'a usize, VarPos<'a>)>;

  fn into_iter(self) -> Self::IntoIter {
    self
      .entries
      .iter()
      .map(|(pos, state, dis)| (pos, VarPos { dis, state: *state }))
  }
}

#[derive(Clone, Debug, Default)]
pub struct VarPosMap {
  entries: Vec<Entry>,
}

impl VarPosMap {
  pub(crate) fn with_capacity(capacity: usize) -> Self {
    Self {
      entries: Vec::with_capacity(capacity),
    }
  }

  pub(crate) fn len(&self) -> usize {
    self.entries.len()
  }

  pub(crate) fn get(&self, pos: usize) -> Option<VarPos<'_>> {
    let index = self.entries.partition_point(|(key, _, _)| *key < pos);
    self
      .entries
      .get(index)
      .filter(|(key, _, _)| *key == pos)
      .map(|(_, state, dis)| VarPos { dis, state: *state })
  }

  pub(crate) fn contains_key(&self, pos: usize) -> bool {
    self.entries.binary_search_by_key(&pos, |(key, _, _)| *key).is_ok()
  }

  pub(crate) fn insert(&mut self, pos: usize, state: AsciiChar, dis: &[f64]) {
    self.insert_with(pos, state, dis.len(), |row| row.copy_from_slice(dis));
  }

  pub(crate) fn insert_with(&mut self, pos: usize, state: AsciiChar, n_states: usize, fill: impl FnOnce(&mut [f64])) {
    let mut dis: SiteProfile = smallvec![0.0; n_states];
    fill(&mut dis);
    let entry = (pos, state, dis);
    if self.entries.last().is_none_or(|(last, _, _)| *last < pos) {
      self.entries.push(entry);
      return;
    }
    match self.entries.binary_search_by_key(&pos, |(key, _, _)| *key) {
      Ok(index) => self.entries[index] = entry,
      Err(index) => self.entries.insert(index, entry),
    }
  }

  pub(crate) fn iter(&self) -> impl Iterator<Item = (&usize, VarPos<'_>)> {
    self.into_iter()
  }

  pub(crate) fn keys(&self) -> impl Iterator<Item = &usize> {
    self.entries.iter().map(|(pos, _, _)| pos)
  }
}

type Entry = (usize, AsciiChar, SiteProfile);
