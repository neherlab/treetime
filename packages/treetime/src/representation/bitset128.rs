use auto_ops::{impl_op_ex, impl_op_ex_commutative};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::borrow::Borrow;

#[allow(variant_size_differences)]
#[derive(Clone, Debug)]
pub enum Bitset128Status {
  Empty,
  Unambiguous(char),
  Ambiguous(BitSet128),
}

#[must_use]
#[derive(Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct BitSet128 {
  bits: u128,
}

impl BitSet128 {
  pub fn new() -> Self {
    Self { bits: 0 }
  }

  pub fn from_char(c: char) -> Self {
    Self { bits: 1 << (c as u32) }
  }

  pub fn from_slice(chars: &[char]) -> Self {
    chars.iter().copied().collect()
  }

  pub fn is_empty(&self) -> bool {
    self.bits == 0
  }

  pub fn len(&self) -> usize {
    self.bits.count_ones() as usize
  }

  pub fn clear(&mut self) {
    self.bits = 0;
  }

  pub fn contains(&self, c: char) -> bool {
    (self.bits & (1 << (c as u32))) != 0
  }

  pub fn insert(&mut self, c: char) {
    let mask = 1 << (c as u32);
    self.bits |= mask;
  }

  pub fn remove(&mut self, c: char) {
    let mask = 1 << (c as u32);
    self.bits &= !mask;
  }

  pub fn union(&self, other: &Self) -> Self {
    Self {
      bits: self.bits | other.bits,
    }
  }

  pub fn intersection(&self, other: &Self) -> Self {
    Self {
      bits: self.bits & other.bits,
    }
  }

  pub fn difference(&self, other: &Self) -> Self {
    Self {
      bits: self.bits & !other.bits,
    }
  }

  pub fn symmetric_difference(&self, other: &Self) -> Self {
    Self {
      bits: self.bits ^ other.bits,
    }
  }

  pub fn from_union<I>(sets: I) -> Self
  where
    I: IntoIterator,
    I::Item: Borrow<Self>,
  {
    let bits = sets.into_iter().fold(0, |acc, set| acc | set.borrow().bits);
    Self { bits }
  }

  pub fn from_intersection<I>(sets: I) -> Self
  where
    I: IntoIterator,
    I::Item: Borrow<Self>,
  {
    sets
      .into_iter()
      .map(|set| set.borrow().bits)
      .reduce(|acc, set| acc & set)
      .map_or_else(Self::new, |bits| Self { bits })
  }

  pub fn is_disjoint(&self, other: &Self) -> bool {
    (self.bits & other.bits) == 0
  }

  pub fn is_subset(&self, other: &Self) -> bool {
    (self.bits & other.bits) == self.bits
  }

  pub fn is_superset(&self, other: &Self) -> bool {
    other.is_subset(self)
  }

  pub fn iter(&self) -> impl Iterator<Item = char> + '_ {
    (0..128)
      .filter(|&i| (self.bits & (1 << i)) != 0)
      .map(|i| char::from_u32(i).unwrap())
  }

  pub fn chars(&self) -> impl Iterator<Item = char> + '_ {
    self.iter()
  }

  pub fn get(&self) -> Bitset128Status {
    match self.bits.count_ones() {
      0 => Bitset128Status::Empty,
      1 => Bitset128Status::Unambiguous(self.get_one()),
      _ => Bitset128Status::Ambiguous(*self),
    }
  }

  pub fn first(&self) -> Option<char> {
    (!self.is_empty()).then_some(char::from_u32(self.bits.trailing_zeros()).unwrap())
  }

  pub fn last(&self) -> Option<char> {
    (!self.is_empty()).then_some(char::from_u32(127 - self.bits.leading_zeros()).unwrap())
  }

  pub fn get_one_maybe(&self) -> Option<char> {
    self.first()
  }

  pub fn get_one(&self) -> char {
    self.get_one_maybe().expect("BitSet128 is empty")
  }

  pub fn get_one_exactly(&self) -> char {
    assert_eq!(self.len(), 1, "expected exactly one element");
    self.get_one()
  }

  pub fn to_vec(&self) -> Vec<char> {
    self.iter().collect()
  }

  pub fn from_vec(chars: Vec<char>) -> Self {
    Self::from_iter(chars)
  }
}

impl<T: Borrow<char>> Extend<T> for BitSet128 {
  fn extend<I>(&mut self, iter: I)
  where
    I: IntoIterator<Item = T>,
  {
    for c in iter {
      self.insert(*c.borrow());
    }
  }
}

impl<T: Borrow<char>> FromIterator<T> for BitSet128 {
  fn from_iter<I>(iter: I) -> Self
  where
    I: IntoIterator<Item = T>,
  {
    let bits = iter.into_iter().fold(0, |acc, c| acc | (1 << (*c.borrow() as u32)));
    Self { bits }
  }
}

impl std::fmt::Display for BitSet128 {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    let chars: String = (0..128)
      .filter(|&i| (self.bits & (1 << i)) != 0)
      .map(|i| char::from_u32(i).unwrap())
      .join(", ");
    write!(f, "{{{chars}}}")
  }
}

impl Serialize for BitSet128 {
  fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
  where
    S: serde::Serializer,
  {
    let chars: Vec<char> = self.to_vec();
    chars.serialize(serializer)
  }
}

impl<'de> Deserialize<'de> for BitSet128 {
  fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
  where
    D: serde::Deserializer<'de>,
  {
    let chars: Vec<char> = Vec::deserialize(deserializer)?;
    Ok(Self::from_vec(chars))
  }
}

impl std::fmt::Debug for BitSet128 {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    write!(f, "{self}")
  }
}

#[macro_export]
macro_rules! bitset128 {
  ($($char:expr),* $(,)?) => {
    {
      let chars = [$($char),*];
      BitSet128::from_slice(&chars)
    }
  };
}

impl_op_ex!(+|a: &BitSet128, b: &BitSet128| -> BitSet128 { a.union(b) });
impl_op_ex!(+=|a: &mut BitSet128, b: &BitSet128| { *a = a.union(b); });

impl_op_ex!(| |a: &BitSet128, b: &BitSet128| -> BitSet128 { a.union(b) });
impl_op_ex!(|= |a: &mut BitSet128, b: &BitSet128| { *a = a.union(b); });

impl_op_ex!(&|a: &BitSet128, b: &BitSet128| -> BitSet128 { a.intersection(b) });
impl_op_ex!(&= |a: &mut BitSet128, b: &BitSet128| { *a = a.intersection(b); });

impl_op_ex!(-|a: &BitSet128, b: &BitSet128| -> BitSet128 { a.difference(b) });
impl_op_ex!(-=|a: &mut BitSet128, b: &BitSet128| { *a = a.difference(b); });

impl_op_ex!(^|a: &BitSet128, b: &BitSet128| -> BitSet128 { a.symmetric_difference(b) });
impl_op_ex!(^=|a: &mut BitSet128, b: &BitSet128| { *a = a.symmetric_difference(b); });

impl_op_ex_commutative!(+|a: &BitSet128, b: char| -> BitSet128 { a.union(&BitSet128::from_char(b)) });
impl_op_ex!(+=|a: &mut BitSet128, b: char| { *a = a.union(&BitSet128::from_char(b)); });

impl_op_ex_commutative!(| |a: &BitSet128, b: char| -> BitSet128 { a.union(&BitSet128::from_char(b)) });
impl_op_ex!(|=|a: &mut BitSet128, b: char| { *a = a.union(&BitSet128::from_char(b)); });

impl_op_ex!(-|a: &BitSet128, b: char| -> BitSet128 { a.difference(&BitSet128::from_char(b)) });
impl_op_ex!(-=|a: &mut BitSet128, b: char| { *a = a.difference(&BitSet128::from_char(b)); });

impl_op_ex_commutative!(&|a: &BitSet128, b: char| -> BitSet128 { a.intersection(&BitSet128::from_char(b)) });
impl_op_ex!(&=|a: &mut BitSet128, b: char| { *a = a.intersection(&BitSet128::from_char(b)); });

impl_op_ex_commutative!(^|a: &BitSet128, b: char| -> BitSet128 { a.symmetric_difference(&BitSet128::from_char(b)) });
impl_op_ex!(^=|a: &mut BitSet128, b: char| { *a = a.symmetric_difference(&BitSet128::from_char(b)); });

#[cfg(test)]
mod tests {
  #![allow(clippy::op_ref)]
  use super::*;
  use pretty_assertions::{assert_eq, assert_ne};
  use rstest::rstest;
  use std::hash::{DefaultHasher, Hash, Hasher};

  #[test]
  fn test_bitset128_new() {
    let actual = BitSet128::new();
    let expected = bitset128! {};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_iter() {
    let actual = BitSet128::from_iter(['a', 'b', 'c']);
    let expected = bitset128! {'a', 'b', 'c'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_slice() {
    let actual = BitSet128::from_slice(&['x', 'y', 'z']);
    let expected = bitset128! {'x', 'y', 'z'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_is_empty() {
    let a = BitSet128::new();
    assert!(a.is_empty());
    let a = bitset128! {'a'};
    assert!(!a.is_empty());
  }

  #[test]
  fn test_bitset128_len() {
    let actual = bitset128! {'a', 'b', 'c'};
    let expected_len = 3;
    assert_eq!(actual.len(), expected_len);

    let actual = BitSet128::new();
    let expected_len = 0;
    assert_eq!(actual.len(), expected_len);
  }

  #[test]
  fn test_bitset128_clear() {
    let mut actual = bitset128! {'a', 'b', 'c'};
    actual.clear();
    let expected = BitSet128::new();
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_insert() {
    let mut actual = bitset128! {};
    actual.insert('a');
    actual.insert('a');
    actual.insert('b');
    let expected = bitset128! {'a', 'b'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_remove() {
    let mut actual = bitset128! {'a', 'b', 'c'};
    actual.remove('b');
    actual.remove('b');
    let expected = bitset128! {'a', 'c'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_union() {
    let a = bitset128! {'a', 'b', 'y'};
    let b = bitset128! {'b', 'z', 'x'};
    let actual = a.union(&b);
    let expected = bitset128! {'a', 'b', 'x', 'y', 'z'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_union_with_empty() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {};
    let actual = a.union(&b);
    let expected = bitset128! {'a', 'b', 'c'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_intersection() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {'b', 'c', 'x'};
    let actual = a.intersection(&b);
    let expected = bitset128! {'b', 'c'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_intersection_with_empty() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {};
    let actual = a.intersection(&b);
    let expected = bitset128! {};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_union() {
    let a = bitset128! {'a', 'b', 'x', 'z'};
    let b = bitset128! {'y', 'x', 'a', 'z'};
    let c = bitset128! {'p', 'q', 'y', 'a'};
    let actual = BitSet128::from_union([a, b, c]);
    let expected = bitset128! {'a', 'b', 'p', 'q', 'x', 'y', 'z'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_union_with_empty() {
    let a = bitset128! {'a', 'b', 'x'};
    let b = bitset128! {};
    let c = bitset128! {'p', 'q', 'y'};
    let actual = BitSet128::from_union([a, b, c]);
    let expected = bitset128! {'a', 'b', 'p', 'q', 'x', 'y'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_union_of_all_empty() {
    let a = bitset128! {};
    let b = bitset128! {};
    let c = bitset128! {};
    let actual = BitSet128::from_union([a, b, c]);
    let expected = bitset128! {};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_intersection() {
    let a = bitset128! {'a', 'b', 'x', 'z'};
    let b = bitset128! {'y', 'x', 'a', 'z'};
    let c = bitset128! {'x', 'q', 'y', 'a'};
    let actual = BitSet128::from_intersection([a, b, c]);
    let expected = bitset128! {'a', 'x'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_intersection_with_empty() {
    let a = bitset128! {'a', 'b', 'x'};
    let b = bitset128! {};
    let c = bitset128! {'p', 'q', 'y'};
    let actual = BitSet128::from_intersection([a, b, c]);
    let expected = bitset128! {};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_intersection_of_all_empty() {
    let a = bitset128! {};
    let b = bitset128! {};
    let c = bitset128! {};
    let actual = BitSet128::from_intersection([a, b, c]);
    let expected = bitset128! {};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_difference() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {'b', 'c', 'x'};
    let actual = a.difference(&b);
    let expected = bitset128! {'a'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_is_disjoint() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {'x', 'y', 'z'};
    assert!(a.is_disjoint(&b));
  }

  #[test]
  fn test_bitset128_is_subset() {
    let a = bitset128! {'a', 'b'};
    let b = bitset128! {'a', 'b', 'c'};
    assert!(a.is_subset(&b));
  }

  #[test]
  fn test_bitset128_is_superset() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {'a', 'b'};
    assert!(a.is_superset(&b));
  }

  #[test]
  fn test_bitset128_display() {
    let a = bitset128! {'a', 'b', 'c'};
    let actual = a.to_string();
    let expected = "{a, b, c}";
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_debug() {
    let a = bitset128! {'a', 'b', 'c'};
    let actual = format!("{a:?}");
    let expected = "{a, b, c}";
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_eq() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {'a', 'b', 'c'};
    let c = bitset128! {'x', 'y', 'z'};
    assert_eq!(a, b);
    assert_ne!(a, c);
  }

  #[test]
  fn test_bitset128_hash() {
    fn calculate_hash<T: Hash>(t: &T) -> u64 {
      let mut hasher = DefaultHasher::new();
      t.hash(&mut hasher);
      hasher.finish()
    }

    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {'a', 'b', 'c'};
    let c = bitset128! {'x', 'y', 'z'};
    assert_eq!(calculate_hash(&a), calculate_hash(&b));
    assert_ne!(calculate_hash(&a), calculate_hash(&c));
  }

  #[test]
  fn test_bitset128_get_empty() {
    let set = BitSet128::new();
    assert!(matches!(set.get(), Bitset128Status::Empty));
  }

  #[test]
  fn test_bitset128_get_unambiguous() {
    let set = BitSet128::from_char('A');
    assert!(matches!(set.get(), Bitset128Status::Unambiguous('A')));
  }

  #[test]
  fn test_bitset128_get_ambiguous() {
    let set = BitSet128::from_iter(vec!['A', 'C']);
    assert!(matches!(set.get(), Bitset128Status::Ambiguous(_)));
    if let Bitset128Status::Ambiguous(actual) = set.get() {
      let expected = bitset128! {'A', 'C'};
      assert_eq!(actual, expected);
    }
  }

  #[test]
  fn test_bitset128_get_one() {
    let set = BitSet128::from_iter(['T', 'A']);
    let actual = set.get_one();
    let expected = 'A';
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_get_one_exactly() {
    let set = BitSet128::from_iter(['T']);
    let actual = set.get_one();
    let expected = 'T';
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q', 'x'})]
  #[case(bitset128!{'z'},            bitset128!{'z'},            bitset128!{'z'})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e', 'f'},  bitset128!{'a', 'b', 'c', 'd', 'e', 'f'})]
  #[trace]
  fn test_bitset128_add(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a + b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q', 'x'})]
  #[case(bitset128!{'z'},            bitset128!{'z'},            bitset128!{'z'})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e', 'f'},  bitset128!{'a', 'b', 'c', 'd', 'e', 'f'})]
  #[trace]
  fn test_bitset128_add_ref_right(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a + &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q', 'x'})]
  #[case(bitset128!{'z'},            bitset128!{'z'},            bitset128!{'z'})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e', 'f'},  bitset128!{'a', 'b', 'c', 'd', 'e', 'f'})]
  #[trace]
  fn test_bitset128_add_ref_left(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a + b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q', 'x'})]
  #[case(bitset128!{'z'},            bitset128!{'z'},            bitset128!{'z'})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e', 'f'},  bitset128!{'a', 'b', 'c', 'd', 'e', 'f'})]
  #[trace]
  fn test_bitset128_add_ref_both(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a + &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q', 'x'})]
  #[case(bitset128!{'z'},            bitset128!{'z'},            bitset128!{'z'})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e', 'f'},  bitset128!{'a', 'b', 'c', 'd', 'e', 'f'})]
  #[trace]
  fn test_bitset128_add_assign(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a += b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q', 'x'})]
  #[case(bitset128!{'z'},            bitset128!{'z'},            bitset128!{'z'})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e', 'f'},  bitset128!{'a', 'b', 'c', 'd', 'e', 'f'})]
  #[trace]
  fn test_bitset128_add_assign_ref(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a += &b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            'b',                        bitset128!{'a', 'b'})]
  #[case(bitset128!{'x', 'y'},       'z',                        bitset128!{'x', 'y', 'z'})]
  #[case(bitset128!{'m', 'n'},       'm',                        bitset128!{'m', 'n'})]
  #[case(bitset128!{},               'a',                        bitset128!{'a'})]
  #[trace]
  fn test_bitset128_add_char_to_set(#[case] a: BitSet128, #[case] b: char, #[case] expected: BitSet128) {
    let actual = a + b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case('a',                        bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case('x',                        bitset128!{'y', 'z'},       bitset128!{'x', 'y', 'z'})]
  #[case('m',                        bitset128!{'m', 'n'},       bitset128!{'m', 'n'})]
  #[case('a',                        bitset128!{},               bitset128!{'a'})]
  #[trace]
  fn test_bitset128_add_set_to_char(#[case] a: char, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = b + a;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'x', 'q'})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'z'},            bitset128!{'x', 'y', 'z'})]
  #[trace]
  fn test_bitset128_or(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a | b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'x', 'q'})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'z'},            bitset128!{'x', 'y', 'z'})]
  #[trace]
  fn test_bitset128_or_ref_right(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a | &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'x', 'q'})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'z'},            bitset128!{'x', 'y', 'z'})]
  #[trace]
  fn test_bitset128_or_ref_left(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a | b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'x', 'q'})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'z'},            bitset128!{'x', 'y', 'z'})]
  #[trace]
  fn test_bitset128_or_ref_both(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a | &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'x', 'q'})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'z'},            bitset128!{'x', 'y', 'z'})]
  #[trace]
  fn test_bitset128_or_assign(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a |= b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'x', 'q'})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'z'},            bitset128!{'x', 'y', 'z'})]
  #[trace]
  fn test_bitset128_or_assign_ref(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a |= &b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            'b',                        bitset128!{'a', 'b'})]
  #[case(bitset128!{'x', 'y'},       'z',                        bitset128!{'x', 'y', 'z'})]
  #[case(bitset128!{'m', 'n'},       'm',                        bitset128!{'m', 'n'})]
  #[case(bitset128!{},               'a',                        bitset128!{'a'})]
  #[trace]
  fn test_bitset128_or_char_to_set(#[case] a: BitSet128, #[case] b: char, #[case] expected: BitSet128) {
    let actual = a | b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case('a',                        bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case('x',                        bitset128!{'y', 'z'},       bitset128!{'x', 'y', 'z'})]
  #[case('m',                        bitset128!{'m', 'n'},       bitset128!{'m', 'n'})]
  #[case('a',                        bitset128!{},               bitset128!{'a'})]
  #[trace]
  fn test_bitset128_or_set_to_char(#[case] a: char, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = b | a;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{})]
  #[case(bitset128!{'a', 'b'},       bitset128!{'b'},            bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'x'})]
  #[trace]
  fn test_bitset128_sub(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a - b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{})]
  #[case(bitset128!{'a', 'b'},       bitset128!{'b'},            bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'x'})]
  #[trace]
  fn test_bitset128_sub_ref_right(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a - &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{})]
  #[case(bitset128!{'a', 'b'},       bitset128!{'b'},            bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'x'})]
  #[trace]
  fn test_bitset128_sub_ref_left(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a - b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{})]
  #[case(bitset128!{'a', 'b'},       bitset128!{'b'},            bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'x'})]
  #[trace]
  fn test_bitset128_sub_ref_both(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a - &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{})]
  #[case(bitset128!{'a', 'b'},       bitset128!{'b'},            bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'x'})]
  #[trace]
  fn test_bitset128_sub_assign(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a -= b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{})]
  #[case(bitset128!{'a', 'b'},       bitset128!{'b'},            bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e'},       bitset128!{'a', 'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'x'})]
  #[trace]
  fn test_bitset128_sub_assign_ref(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a -= &b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case(bitset128!{},               'b',                        bitset128!{})]
  #[case(bitset128!{'a'},            'b',                        bitset128!{'a'})]
  #[case(bitset128!{'x', 'y'},       'y',                        bitset128!{'x'})]
  #[case(bitset128!{'m', 'n'},       'm',                        bitset128!{'n'})]
  #[case(bitset128!{'a'},            'a',                        bitset128!{})]
  #[trace]
  fn test_bitset128_sub_char_from_set(#[case] a: BitSet128, #[case] b: char, #[case] expected: BitSet128) {
    let actual = a - b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'x'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'a', 'd', 'e'},  bitset128!{'a'})]
  #[case(bitset128!{'m', 'n', 'o'},  bitset128!{'n', 'o', 'p'},  bitset128!{'n', 'o'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'y'})]
  #[case(bitset128!{'g', 'h', 'i'},  bitset128!{'a', 'b', 'c'},  bitset128!{})]
  #[trace]
  fn test_bitset128_and(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a & b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'x'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'a', 'd', 'e'},  bitset128!{'a'})]
  #[case(bitset128!{'m', 'n', 'o'},  bitset128!{'n', 'o', 'p'},  bitset128!{'n', 'o'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'y'})]
  #[case(bitset128!{'g', 'h', 'i'},  bitset128!{'a', 'b', 'c'},  bitset128!{})]
  #[trace]
  fn test_bitset128_and_ref_right(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a & &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'x'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'a', 'd', 'e'},  bitset128!{'a'})]
  #[case(bitset128!{'m', 'n', 'o'},  bitset128!{'n', 'o', 'p'},  bitset128!{'n', 'o'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'y'})]
  #[case(bitset128!{'g', 'h', 'i'},  bitset128!{'a', 'b', 'c'},  bitset128!{})]
  #[trace]
  fn test_bitset128_and_ref_left(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a & b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'x'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'a', 'd', 'e'},  bitset128!{'a'})]
  #[case(bitset128!{'m', 'n', 'o'},  bitset128!{'n', 'o', 'p'},  bitset128!{'n', 'o'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'y'})]
  #[case(bitset128!{'g', 'h', 'i'},  bitset128!{'a', 'b', 'c'},  bitset128!{})]
  #[trace]
  fn test_bitset128_and_ref_both(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a & &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'x'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'a', 'd', 'e'},  bitset128!{'a'})]
  #[case(bitset128!{'m', 'n', 'o'},  bitset128!{'n', 'o', 'p'},  bitset128!{'n', 'o'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'y'})]
  #[case(bitset128!{'g', 'h', 'i'},  bitset128!{'a', 'b', 'c'},  bitset128!{})]
  #[trace]
  fn test_bitset128_and_assign(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a &= b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{})]
  #[case(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'b', 'c'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'x'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'a', 'd', 'e'},  bitset128!{'a'})]
  #[case(bitset128!{'m', 'n', 'o'},  bitset128!{'n', 'o', 'p'},  bitset128!{'n', 'o'})]
  #[case(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'y'})]
  #[case(bitset128!{'g', 'h', 'i'},  bitset128!{'a', 'b', 'c'},  bitset128!{})]
  #[trace]
  fn test_bitset128_and_assign_ref(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a &= &b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case(bitset128!{'a', 'b'},       'a',                        bitset128!{'a'})]
  #[case(bitset128!{'x', 'y'},       'z',                        bitset128!{})]
  #[case(bitset128!{'m', 'n', 'o'},  'o',                        bitset128!{'o'})]
  #[trace]
  fn test_bitset128_and_char_to_set(#[case] a: BitSet128, #[case] b: char, #[case] expected: BitSet128) {
    let actual = a & b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case('a',                        bitset128!{'b', 'a'},       bitset128!{'a'})]
  #[case('z',                        bitset128!{'x', 'y'},       bitset128!{})]
  #[case('o',                        bitset128!{'m', 'n', 'o'},  bitset128!{'o'})]
  #[trace]
  fn test_bitset128_and_set_to_char(#[case] a: char, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = b & a;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q'})]
  #[case(bitset128!{'a', 'b'},       bitset128!{'a', 'b'},       bitset128!{})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'n', 'o'},       bitset128!{'m', 'o'})]
  #[case(bitset128!{'x', 'y', 'z'},  bitset128!{'y', 'z', 'a'},  bitset128!{'x', 'a'})]
  #[case(bitset128!{'a', 'c', 'e'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a', 'b', 'd', 'e'})]
  #[trace]
  fn test_bitset128_xor(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a ^ b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q'})]
  #[case(bitset128!{'a', 'b'},       bitset128!{'a', 'b'},       bitset128!{})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'n', 'o'},       bitset128!{'m', 'o'})]
  #[case(bitset128!{'x', 'y', 'z'},  bitset128!{'y', 'z', 'a'},  bitset128!{'x', 'a'})]
  #[case(bitset128!{'a', 'c', 'e'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a', 'b', 'd', 'e'})]
  #[trace]
  fn test_bitset128_xor_ref_right(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a ^ &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q'})]
  #[case(bitset128!{'a', 'b'},       bitset128!{'a', 'b'},       bitset128!{})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'n', 'o'},       bitset128!{'m', 'o'})]
  #[case(bitset128!{'x', 'y', 'z'},  bitset128!{'y', 'z', 'a'},  bitset128!{'x', 'a'})]
  #[case(bitset128!{'a', 'c', 'e'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a', 'b', 'd', 'e'})]
  #[trace]
  fn test_bitset128_xor_ref_left(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a ^ b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q'})]
  #[case(bitset128!{'a', 'b'},       bitset128!{'a', 'b'},       bitset128!{})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'n', 'o'},       bitset128!{'m', 'o'})]
  #[case(bitset128!{'x', 'y', 'z'},  bitset128!{'y', 'z', 'a'},  bitset128!{'x', 'a'})]
  #[case(bitset128!{'a', 'c', 'e'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a', 'b', 'd', 'e'})]
  #[trace]
  fn test_bitset128_xor_ref_both(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a ^ &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q'})]
  #[case(bitset128!{'a', 'b'},       bitset128!{'a', 'b'},       bitset128!{})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'n', 'o'},       bitset128!{'m', 'o'})]
  #[case(bitset128!{'x', 'y', 'z'},  bitset128!{'y', 'z', 'a'},  bitset128!{'x', 'a'})]
  #[case(bitset128!{'a', 'c', 'e'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a', 'b', 'd', 'e'})]
  #[trace]
  fn test_bitset128_xor_assign(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a ^= b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a'})]
  #[case(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q'})]
  #[case(bitset128!{'a', 'b'},       bitset128!{'a', 'b'},       bitset128!{})]
  #[case(bitset128!{'m', 'n'},       bitset128!{'n', 'o'},       bitset128!{'m', 'o'})]
  #[case(bitset128!{'x', 'y', 'z'},  bitset128!{'y', 'z', 'a'},  bitset128!{'x', 'a'})]
  #[case(bitset128!{'a', 'c', 'e'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a', 'b', 'd', 'e'})]
  #[trace]
  fn test_bitset128_xor_assign_ref(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a ^= &b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case(bitset128!{'a'},            'b',                        bitset128!{'a', 'b'})]
  #[case(bitset128!{'x', 'y'},       'x',                        bitset128!{'y'})]
  #[case(bitset128!{'m', 'n'},       'm',                        bitset128!{'n'})]
  #[case(bitset128!{'c', 'd'},       'd',                        bitset128!{'c'})]
  #[trace]
  fn test_bitset128_xor_char_to_set(#[case] a: BitSet128, #[case] b: char, #[case] expected: BitSet128) {
    let actual = a ^ b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case('a',                        bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case('x',                        bitset128!{'x', 'y'},       bitset128!{'y'})]
  #[case('m',                        bitset128!{'m', 'n'},       bitset128!{'n'})]
  #[case('d',                        bitset128!{'c', 'd'},       bitset128!{'c'})]
  #[trace]
  fn test_bitset128_xor_set_to_char(#[case] a: char, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = b ^ a;
    assert_eq!(actual, expected);
  }
}
