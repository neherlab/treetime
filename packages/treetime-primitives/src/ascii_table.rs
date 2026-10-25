use crate::seq_char::AsciiChar;
use std::array;
use std::ops::{Index, IndexMut};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AsciiTable<T>([T; AsciiChar::COUNT]);

impl<T> AsciiTable<T> {
  pub fn from_fn(mut entry: impl FnMut(AsciiChar) -> T) -> Self {
    let mut byte = 0_u8;
    Self(array::from_fn(|_| {
      let c = AsciiChar::from_byte_unchecked(byte);
      byte += 1;
      entry(c)
    }))
  }

  pub fn iter(&self) -> impl Iterator<Item = (AsciiChar, &T)> {
    (0_u8..)
      .zip(&self.0)
      .map(|(byte, value)| (AsciiChar::from_byte_unchecked(byte), value))
  }
}

impl<T: Default> Default for AsciiTable<T> {
  fn default() -> Self {
    Self::from_fn(|_| T::default())
  }
}

impl<T> Index<AsciiChar> for AsciiTable<T> {
  type Output = T;

  fn index(&self, c: AsciiChar) -> &T {
    &self.0[usize::from(c)]
  }
}

impl<T> IndexMut<AsciiChar> for AsciiTable<T> {
  fn index_mut(&mut self, c: AsciiChar) -> &mut T {
    &mut self.0[usize::from(c)]
  }
}
