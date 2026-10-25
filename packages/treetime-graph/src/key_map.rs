use crate::edge::GraphEdgeKey;
use crate::node::GraphNodeKey;
use rayon::iter::{FromParallelIterator, IntoParallelIterator, ParallelIterator};
use std::borrow::Borrow;
use std::fmt::{self, Debug, Display};
use std::iter::{Enumerate, FusedIterator};
use std::marker::PhantomData;
use std::ops::Index;
use std::{slice, vec};
use treetime_utils::testing::map_like::MapLike;

#[macro_export]
macro_rules! keymap {
  ($($key:expr => $value:expr),* $(,)?) => {
    <$crate::key_map::KeyMap<_, _> as ::core::iter::FromIterator<_>>::from_iter([$(($key, $value)),*])
  };
}

pub type NodeMap<V> = KeyMap<GraphNodeKey, V>;

pub type EdgeMap<V> = KeyMap<GraphEdgeKey, V>;

#[derive(Clone)]
pub struct KeyMap<K, V> {
  slots: Vec<Option<V>>,
  len: usize,
  key: PhantomData<fn() -> K>,
}

impl<K, V> KeyMap<K, V> {
  pub const fn new() -> Self {
    Self {
      slots: Vec::new(),
      len: 0,
      key: PhantomData,
    }
  }

  pub const fn len(&self) -> usize {
    self.len
  }

  pub const fn is_empty(&self) -> bool {
    self.len == 0
  }

  pub fn values(&self) -> impl DoubleEndedIterator<Item = &V> {
    self.slots.iter().flatten()
  }

  pub fn values_mut(&mut self) -> impl DoubleEndedIterator<Item = &mut V> {
    self.slots.iter_mut().flatten()
  }

  pub fn into_values(self) -> impl DoubleEndedIterator<Item = V> {
    self.slots.into_iter().flatten()
  }

  pub(crate) fn from_slots(slots: Vec<Option<V>>) -> Self {
    Self {
      len: slots.iter().flatten().count(),
      slots,
      key: PhantomData,
    }
  }
}

impl<K: GraphKey, V> KeyMap<K, V> {
  pub fn get(&self, key: impl Borrow<K>) -> Option<&V> {
    self.slots.get(key.borrow().slot())?.as_ref()
  }

  pub fn get_mut(&mut self, key: impl Borrow<K>) -> Option<&mut V> {
    self.slots.get_mut(key.borrow().slot())?.as_mut()
  }

  pub fn contains_key(&self, key: impl Borrow<K>) -> bool {
    self.get(key).is_some()
  }

  pub fn insert(&mut self, key: K, value: V) -> Option<V> {
    let previous = self.slot_mut(key).replace(value);
    if previous.is_none() {
      self.len += 1;
    }
    previous
  }

  pub fn remove(&mut self, key: impl Borrow<K>) -> Option<V> {
    let removed = self.slots.get_mut(key.borrow().slot())?.take();
    if removed.is_some() {
      self.len -= 1;
    }
    removed
  }

  pub fn entry(&mut self, key: K) -> Entry<'_, K, V> {
    Entry { map: self, key }
  }

  pub fn retain(&mut self, mut keep: impl FnMut(K, &mut V) -> bool) {
    for (index, slot) in self.slots.iter_mut().enumerate() {
      if let Some(value) = slot
        && !keep(K::from_slot(index, KeyToken(())), value)
      {
        *slot = None;
        self.len -= 1;
      }
    }
  }

  pub fn iter(&self) -> Iter<'_, K, V> {
    Iter {
      slots: self.slots.iter().enumerate(),
      remaining: self.len,
      key: PhantomData,
    }
  }

  pub fn iter_mut(&mut self) -> IterMut<'_, K, V> {
    IterMut {
      slots: self.slots.iter_mut().enumerate(),
      remaining: self.len,
      key: PhantomData,
    }
  }

  pub fn keys(&self) -> impl DoubleEndedIterator<Item = K> + ExactSizeIterator {
    self.iter().map(|(key, _)| key)
  }

  pub fn into_keys(self) -> impl ExactSizeIterator<Item = K> {
    self.into_iter().map(|(key, _)| key)
  }

  fn get_or_insert_with(&mut self, key: K, value: impl FnOnce() -> V) -> &mut V {
    if self.get(key).is_none() {
      self.len += 1;
    }
    self.slot_mut(key).get_or_insert_with(value)
  }

  fn slot_mut(&mut self, key: K) -> &mut Option<V> {
    let slot = key.slot();
    if slot >= self.slots.len() {
      self.slots.resize_with(slot + 1, || None);
    }
    &mut self.slots[slot]
  }
}

pub trait GraphKey: Copy + Ord + Display + Debug {
  const KIND: &'static str;

  fn slot(self) -> usize;

  fn from_slot(slot: usize, token: KeyToken) -> Self;
}

pub struct KeyToken(());

impl GraphKey for GraphNodeKey {
  const KIND: &'static str = "Node";

  fn slot(self) -> usize {
    self.0
  }

  fn from_slot(slot: usize, _token: KeyToken) -> Self {
    Self(slot)
  }
}

impl GraphKey for GraphEdgeKey {
  const KIND: &'static str = "Edge";

  fn slot(self) -> usize {
    self.0
  }

  fn from_slot(slot: usize, _token: KeyToken) -> Self {
    Self(slot)
  }
}

pub struct Entry<'a, K, V> {
  map: &'a mut KeyMap<K, V>,
  key: K,
}

impl<'a, K: GraphKey, V> Entry<'a, K, V> {
  pub fn or_insert_with(self, value: impl FnOnce() -> V) -> &'a mut V {
    self.map.get_or_insert_with(self.key, value)
  }

  pub fn or_insert(self, value: V) -> &'a mut V {
    self.or_insert_with(|| value)
  }

  pub fn or_default(self) -> &'a mut V
  where
    V: Default,
  {
    self.or_insert_with(V::default)
  }

  #[must_use]
  pub fn and_modify(self, modify: impl FnOnce(&mut V)) -> Self {
    if let Some(value) = self.map.get_mut(self.key) {
      modify(value);
    }
    self
  }
}

pub struct Iter<'a, K, V> {
  slots: Enumerate<slice::Iter<'a, Option<V>>>,
  remaining: usize,
  key: PhantomData<fn() -> K>,
}

impl<'a, K: GraphKey, V> Iterator for Iter<'a, K, V> {
  type Item = (K, &'a V);

  fn next(&mut self) -> Option<Self::Item> {
    let entry = self
      .slots
      .find_map(|(index, slot)| Some((K::from_slot(index, KeyToken(())), slot.as_ref()?)))?;
    self.remaining -= 1;
    Some(entry)
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    (self.remaining, Some(self.remaining))
  }
}

impl<K: GraphKey, V> DoubleEndedIterator for Iter<'_, K, V> {
  fn next_back(&mut self) -> Option<Self::Item> {
    let entry = self
      .slots
      .rfind(|(_, slot)| slot.is_some())
      .and_then(|(index, slot)| Some((K::from_slot(index, KeyToken(())), slot.as_ref()?)))?;
    self.remaining -= 1;
    Some(entry)
  }
}

impl<K: GraphKey, V> ExactSizeIterator for Iter<'_, K, V> {}

impl<K: GraphKey, V> FusedIterator for Iter<'_, K, V> {}

pub struct IterMut<'a, K, V> {
  slots: Enumerate<slice::IterMut<'a, Option<V>>>,
  remaining: usize,
  key: PhantomData<fn() -> K>,
}

impl<'a, K: GraphKey, V> Iterator for IterMut<'a, K, V> {
  type Item = (K, &'a mut V);

  fn next(&mut self) -> Option<Self::Item> {
    let entry = self
      .slots
      .find_map(|(index, slot)| Some((K::from_slot(index, KeyToken(())), slot.as_mut()?)))?;
    self.remaining -= 1;
    Some(entry)
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    (self.remaining, Some(self.remaining))
  }
}

impl<K: GraphKey, V> ExactSizeIterator for IterMut<'_, K, V> {}

impl<K: GraphKey, V> FusedIterator for IterMut<'_, K, V> {}

pub struct IntoIter<K, V> {
  slots: Enumerate<vec::IntoIter<Option<V>>>,
  remaining: usize,
  key: PhantomData<fn() -> K>,
}

impl<K: GraphKey, V> Iterator for IntoIter<K, V> {
  type Item = (K, V);

  fn next(&mut self) -> Option<Self::Item> {
    let entry = self
      .slots
      .find_map(|(index, slot)| Some((K::from_slot(index, KeyToken(())), slot?)))?;
    self.remaining -= 1;
    Some(entry)
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    (self.remaining, Some(self.remaining))
  }
}

impl<K: GraphKey, V> ExactSizeIterator for IntoIter<K, V> {}

impl<K: GraphKey, V> FusedIterator for IntoIter<K, V> {}

impl<'a, K: GraphKey, V> IntoIterator for &'a KeyMap<K, V> {
  type Item = (K, &'a V);
  type IntoIter = Iter<'a, K, V>;

  fn into_iter(self) -> Self::IntoIter {
    self.iter()
  }
}

impl<'a, K: GraphKey, V> IntoIterator for &'a mut KeyMap<K, V> {
  type Item = (K, &'a mut V);
  type IntoIter = IterMut<'a, K, V>;

  fn into_iter(self) -> Self::IntoIter {
    self.iter_mut()
  }
}

impl<K: GraphKey, V> IntoIterator for KeyMap<K, V> {
  type Item = (K, V);
  type IntoIter = IntoIter<K, V>;

  fn into_iter(self) -> Self::IntoIter {
    IntoIter {
      slots: self.slots.into_iter().enumerate(),
      remaining: self.len,
      key: PhantomData,
    }
  }
}

impl<K: GraphKey, V> FromIterator<(K, V)> for KeyMap<K, V> {
  fn from_iter<I: IntoIterator<Item = (K, V)>>(entries: I) -> Self {
    let mut map = Self::new();
    map.extend(entries);
    map
  }
}

impl<K: GraphKey, V, const N: usize> From<[(K, V); N]> for KeyMap<K, V> {
  fn from(entries: [(K, V); N]) -> Self {
    entries.into_iter().collect()
  }
}

impl<K: GraphKey + Send, V: Send> FromParallelIterator<(K, V)> for KeyMap<K, V> {
  fn from_par_iter<I: IntoParallelIterator<Item = (K, V)>>(entries: I) -> Self {
    entries.into_par_iter().collect::<Vec<_>>().into_iter().collect()
  }
}

impl<K: GraphKey, V> Extend<(K, V)> for KeyMap<K, V> {
  fn extend<I: IntoIterator<Item = (K, V)>>(&mut self, entries: I) {
    for (key, value) in entries {
      self.insert(key, value);
    }
  }
}

impl<K: GraphKey, V> Index<K> for KeyMap<K, V> {
  type Output = V;

  fn index(&self, key: K) -> &V {
    self.get(key).unwrap_or_else(|| missing_entry(K::KIND, key))
  }
}

impl<K: GraphKey, V> Index<&K> for KeyMap<K, V> {
  type Output = V;

  fn index(&self, key: &K) -> &V {
    &self[*key]
  }
}

impl<K: GraphKey, V> MapLike<K, V> for KeyMap<K, V> {
  fn map_get(&self, key: &K) -> Option<&V> {
    self.get(key)
  }

  fn map_len(&self) -> usize {
    self.len
  }

  fn map_iter<'a>(&'a self) -> impl Iterator<Item = (K, &'a V)>
  where
    K: 'a,
    V: 'a,
  {
    self.iter()
  }
}

impl<K, V> Default for KeyMap<K, V> {
  fn default() -> Self {
    Self::new()
  }
}

impl<K: GraphKey, V: PartialEq> PartialEq for KeyMap<K, V> {
  fn eq(&self, other: &Self) -> bool {
    self.len == other.len && self.iter().eq(other.iter())
  }
}

impl<K: GraphKey, V: Eq> Eq for KeyMap<K, V> {}

impl<K: GraphKey, V: Debug> Debug for KeyMap<K, V> {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    f.debug_map().entries(self.iter()).finish()
  }
}

#[expect(
  clippy::panic,
  reason = "a key without an entry is a programming error, as an absent key is for a map index"
)]
fn missing_entry(kind: &str, key: impl Display) -> ! {
  panic!("{kind} {key} has no entry in the map")
}
