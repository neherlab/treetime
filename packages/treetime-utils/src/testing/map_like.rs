use indexmap::IndexMap;
use std::collections::{BTreeMap, HashMap};
use std::hash::Hash;

pub trait MapLike<K, V> {
  fn map_get(&self, key: &K) -> Option<&V>;
  fn map_len(&self) -> usize;
  fn map_iter<'a>(&'a self) -> impl Iterator<Item = (&'a K, &'a V)>
  where
    K: 'a,
    V: 'a;
}

impl<K: Ord, V> MapLike<K, V> for BTreeMap<K, V> {
  fn map_get(&self, key: &K) -> Option<&V> {
    self.get(key)
  }
  fn map_len(&self) -> usize {
    self.len()
  }
  fn map_iter<'a>(&'a self) -> impl Iterator<Item = (&'a K, &'a V)>
  where
    K: 'a,
    V: 'a,
  {
    self.iter()
  }
}

impl<K: Hash + Eq, V> MapLike<K, V> for HashMap<K, V> {
  fn map_get(&self, key: &K) -> Option<&V> {
    self.get(key)
  }
  fn map_len(&self) -> usize {
    self.len()
  }
  fn map_iter<'a>(&'a self) -> impl Iterator<Item = (&'a K, &'a V)>
  where
    K: 'a,
    V: 'a,
  {
    self.iter()
  }
}

impl<K: Hash + Eq, V> MapLike<K, V> for IndexMap<K, V> {
  fn map_get(&self, key: &K) -> Option<&V> {
    self.get(key)
  }
  fn map_len(&self) -> usize {
    self.len()
  }
  fn map_iter<'a>(&'a self) -> impl Iterator<Item = (&'a K, &'a V)>
  where
    K: 'a,
    V: 'a,
  {
    self.iter()
  }
}

impl<K, V, M: MapLike<K, V>> MapLike<K, V> for &M {
  fn map_get(&self, key: &K) -> Option<&V> {
    (*self).map_get(key)
  }
  fn map_len(&self) -> usize {
    (*self).map_len()
  }
  fn map_iter<'a>(&'a self) -> impl Iterator<Item = (&'a K, &'a V)>
  where
    K: 'a,
    V: 'a,
  {
    (*self).map_iter()
  }
}
