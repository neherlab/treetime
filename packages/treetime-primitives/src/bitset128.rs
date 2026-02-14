use crate::seq_char::AsciiChar;
use auto_ops::{impl_op_ex, impl_op_ex_commutative};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::borrow::Borrow;

#[allow(variant_size_differences)]
#[derive(Clone, Debug)]
pub enum Bitset128Status {
  Empty,
  Unambiguous(AsciiChar),
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

  /// Creates a bitset with a single bit set at position `c`.
  ///
  /// # Panics (debug only)
  /// Panics if `c >= 128`.
  pub fn from_char<T: Into<u32>>(c: T) -> Self {
    let c = c.into();
    debug_assert!(c < 128, "BitSet128::from_char requires c < 128, got {c}");
    Self { bits: 1 << c }
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

  /// Returns true if the bit at position `c` is set.
  ///
  /// # Panics (debug only)
  /// Panics if `c >= 128`.
  pub fn contains<T: Into<u32>>(&self, c: T) -> bool {
    let c = c.into();
    debug_assert!(c < 128, "BitSet128::contains requires c < 128, got {c}");
    (self.bits & (1 << c)) != 0
  }

  /// Sets the bit at position `c`.
  ///
  /// # Panics (debug only)
  /// Panics if `c >= 128`.
  pub fn insert<T: Into<u32>>(&mut self, c: T) {
    let c = c.into();
    debug_assert!(c < 128, "BitSet128::insert requires c < 128, got {c}");
    self.bits |= 1 << c;
  }

  /// Clears the bit at position `c`.
  ///
  /// # Panics (debug only)
  /// Panics if `c >= 128`.
  pub fn remove<T: Into<u32>>(&mut self, c: T) {
    let c = c.into();
    debug_assert!(c < 128, "BitSet128::remove requires c < 128, got {c}");
    self.bits &= !(1 << c);
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

  pub fn iter(&self) -> impl Iterator<Item = AsciiChar> + '_ {
    (0..128).filter(|&i| (self.bits & (1 << i)) != 0).map(AsciiChar)
  }

  pub fn chars(&self) -> impl Iterator<Item = AsciiChar> + '_ {
    self.iter()
  }

  pub fn get(&self) -> Bitset128Status {
    match self.bits.count_ones() {
      0 => Bitset128Status::Empty,
      1 => Bitset128Status::Unambiguous(self.get_one()),
      _ => Bitset128Status::Ambiguous(*self),
    }
  }

  pub fn first(&self) -> Option<AsciiChar> {
    (!self.is_empty()).then_some((self.bits.trailing_zeros() as u8).into())
  }

  pub fn last(&self) -> Option<AsciiChar> {
    (!self.is_empty()).then_some(((127 - self.bits.leading_zeros()) as u8).into())
  }

  pub fn get_one_maybe(&self) -> Option<AsciiChar> {
    self.first()
  }

  pub fn get_one(&self) -> AsciiChar {
    self.get_one_maybe().expect("BitSet128 is empty")
  }

  pub fn get_one_exactly(&self) -> AsciiChar {
    assert_eq!(self.len(), 1, "expected exactly one element");
    self.get_one()
  }

  pub fn to_vec(&self) -> Vec<AsciiChar> {
    self.iter().collect()
  }
}

impl<T: Borrow<u32>> Extend<T> for BitSet128 {
  fn extend<I>(&mut self, iter: I)
  where
    I: IntoIterator<Item = T>,
  {
    for c in iter {
      self.insert(*c.borrow());
    }
  }
}

impl FromIterator<AsciiChar> for BitSet128 {
  fn from_iter<I: IntoIterator<Item = AsciiChar>>(iter: I) -> Self {
    let bits = iter.into_iter().fold(0_u128, |acc, c| acc | (1 << c.inner()));
    Self { bits }
  }
}

impl<'a> FromIterator<&'a AsciiChar> for BitSet128 {
  fn from_iter<I: IntoIterator<Item = &'a AsciiChar>>(iter: I) -> Self {
    iter.into_iter().copied().collect()
  }
}

impl FromIterator<u8> for BitSet128 {
  fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
    let bits = iter.into_iter().fold(0_u128, |acc, c| acc | (1 << c));
    Self { bits }
  }
}

impl<'a> FromIterator<&'a u8> for BitSet128 {
  fn from_iter<I: IntoIterator<Item = &'a u8>>(iter: I) -> Self {
    iter.into_iter().copied().collect()
  }
}

impl FromIterator<char> for BitSet128 {
  fn from_iter<I: IntoIterator<Item = char>>(iter: I) -> Self {
    iter.into_iter().map(|c| c as u8).collect()
  }
}

impl<'a> FromIterator<&'a char> for BitSet128 {
  fn from_iter<I: IntoIterator<Item = &'a char>>(iter: I) -> Self {
    iter.into_iter().map(|&c| c as u8).collect()
  }
}

impl From<&[u8]> for BitSet128 {
  fn from(slice: &[u8]) -> Self {
    slice.iter().copied().collect()
  }
}

impl From<Vec<u8>> for BitSet128 {
  fn from(vec: Vec<u8>) -> Self {
    vec.into_iter().collect()
  }
}

impl From<&[char]> for BitSet128 {
  fn from(slice: &[char]) -> Self {
    slice.iter().copied().collect()
  }
}

impl From<Vec<char>> for BitSet128 {
  fn from(vec: Vec<char>) -> Self {
    vec.into_iter().collect()
  }
}

impl From<u8> for BitSet128 {
  fn from(c: u8) -> Self {
    BitSet128 { bits: 1 << c }
  }
}

impl From<char> for BitSet128 {
  fn from(c: char) -> Self {
    BitSet128 { bits: 1 << (c as u8) }
  }
}

impl std::fmt::Display for BitSet128 {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    let chars: String = (0..128)
      .filter(|&i| (self.bits & (1 << i)) != 0)
      .map(|c| c as u8 as char)
      .join(", ");
    write!(f, "{{{chars}}}")
  }
}

impl Serialize for BitSet128 {
  fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
  where
    S: serde::Serializer,
  {
    let chars: Vec<char> = self.iter().map(char::from).collect();
    chars.serialize(serializer)
  }
}

impl<'de> Deserialize<'de> for BitSet128 {
  fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
  where
    D: serde::Deserializer<'de>,
  {
    let chars: Vec<char> = Vec::deserialize(deserializer)?;
    Ok(Self::from_iter(chars))
  }
}

impl std::fmt::Debug for BitSet128 {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    write!(f, "{self}")
  }
}

#[macro_export]
macro_rules! bitset128 {
  () => {
    {
      BitSet128::new()
    }
  };
  ($($char:expr),* $(,)?) => {
    {
      BitSet128::from_iter(&[$($char),*])
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

impl_op_ex_commutative!(+|a: &BitSet128, b: u8| -> BitSet128 { a.union(&BitSet128::from_char(b)) });
impl_op_ex!(+=|a: &mut BitSet128, b: u8| { *a = a.union(&BitSet128::from_char(b)); });

impl_op_ex_commutative!(| |a: &BitSet128, b: u8| -> BitSet128 { a.union(&BitSet128::from_char(b)) });
impl_op_ex!(|=|a: &mut BitSet128, b: u8| { *a = a.union(&BitSet128::from_char(b)); });

impl_op_ex!(-|a: &BitSet128, b: u8| -> BitSet128 { a.difference(&BitSet128::from_char(b)) });
impl_op_ex!(-=|a: &mut BitSet128, b: u8| { *a = a.difference(&BitSet128::from_char(b)); });

impl_op_ex_commutative!(&|a: &BitSet128, b: u8| -> BitSet128 { a.intersection(&BitSet128::from_char(b)) });
impl_op_ex!(&=|a: &mut BitSet128, b: u8| { *a = a.intersection(&BitSet128::from_char(b)); });

impl_op_ex_commutative!(^|a: &BitSet128, b: u8| -> BitSet128 { a.symmetric_difference(&BitSet128::from_char(b)) });
impl_op_ex!(^=|a: &mut BitSet128, b: u8| { *a = a.symmetric_difference(&BitSet128::from_char(b)); });

impl_op_ex_commutative!(+|a: &BitSet128, b: AsciiChar| -> BitSet128 { a.union(&BitSet128::from_char(b)) });
impl_op_ex!(+=|a: &mut BitSet128, b: AsciiChar| { *a = a.union(&BitSet128::from_char(b)); });

impl_op_ex_commutative!(| |a: &BitSet128, b: AsciiChar| -> BitSet128 { a.union(&BitSet128::from_char(b)) });
impl_op_ex!(|=|a: &mut BitSet128, b: AsciiChar| { *a = a.union(&BitSet128::from_char(b)); });

impl_op_ex!(-|a: &BitSet128, b: AsciiChar| -> BitSet128 { a.difference(&BitSet128::from_char(b)) });
impl_op_ex!(-=|a: &mut BitSet128, b: AsciiChar| { *a = a.difference(&BitSet128::from_char(b)); });

impl_op_ex_commutative!(&|a: &BitSet128, b: AsciiChar| -> BitSet128 { a.intersection(&BitSet128::from_char(b)) });
impl_op_ex!(&=|a: &mut BitSet128, b: AsciiChar| { *a = a.intersection(&BitSet128::from_char(b)); });

impl_op_ex_commutative!(^|a: &BitSet128, b: AsciiChar| -> BitSet128 { a.symmetric_difference(&BitSet128::from_char(b)) });
impl_op_ex!(^=|a: &mut BitSet128, b: AsciiChar| { *a = a.symmetric_difference(&BitSet128::from_char(b)); });
