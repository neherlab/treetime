use serde::{Deserialize, Serialize};
use std::fmt::Write as StdFmtWrite;

#[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct AsciiChar(pub u8);

impl AsciiChar {
  pub const fn new(value: u8) -> Self {
    Self(value)
  }

  pub const fn inner(&self) -> u8 {
    self.0
  }
}

impl core::fmt::Display for AsciiChar {
  fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
    f.write_char(self.0 as char)
  }
}

impl core::fmt::Debug for AsciiChar {
  fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
    core::fmt::Display::fmt(self, f)
  }
}

impl From<u8> for AsciiChar {
  fn from(item: u8) -> Self {
    AsciiChar(item)
  }
}

/// # Panics
///
/// Debug builds panic if `item >= 256`.
impl From<u16> for AsciiChar {
  fn from(item: u16) -> Self {
    debug_assert!(item < 256, "AsciiChar::from(u16): value {item} >= 256");
    AsciiChar(item as u8)
  }
}

/// # Panics
///
/// Debug builds panic if `item >= 256`.
impl From<u32> for AsciiChar {
  fn from(item: u32) -> Self {
    debug_assert!(item < 256, "AsciiChar::from(u32): value {item} >= 256");
    AsciiChar(item as u8)
  }
}

/// # Panics
///
/// Debug builds panic if `item >= 256`.
impl From<u64> for AsciiChar {
  fn from(item: u64) -> Self {
    debug_assert!(item < 256, "AsciiChar::from(u64): value {item} >= 256");
    AsciiChar(item as u8)
  }
}

/// # Panics
///
/// Debug builds panic if `item >= 256`.
impl From<usize> for AsciiChar {
  fn from(item: usize) -> Self {
    debug_assert!(item < 256, "AsciiChar::from(usize): value {item} >= 256");
    AsciiChar(item as u8)
  }
}

/// # Panics
///
/// Debug builds panic if `item` is not ASCII (code point >= 128).
impl From<char> for AsciiChar {
  fn from(item: char) -> Self {
    debug_assert!(item.is_ascii(), "AsciiChar::from(char): '{item}' is not ASCII");
    AsciiChar(item as u8)
  }
}

impl From<AsciiChar> for u8 {
  fn from(item: AsciiChar) -> Self {
    item.0
  }
}

impl From<AsciiChar> for u16 {
  fn from(item: AsciiChar) -> Self {
    item.0 as u16
  }
}

impl From<AsciiChar> for u32 {
  fn from(item: AsciiChar) -> Self {
    item.0 as u32
  }
}

impl From<AsciiChar> for u64 {
  fn from(item: AsciiChar) -> Self {
    item.0 as u64
  }
}

impl From<AsciiChar> for usize {
  fn from(item: AsciiChar) -> Self {
    item.0 as usize
  }
}

impl From<AsciiChar> for char {
  fn from(item: AsciiChar) -> Self {
    item.0 as char
  }
}
