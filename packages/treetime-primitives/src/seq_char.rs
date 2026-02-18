use serde::de::{self, Visitor};
use serde::{Deserialize, Deserializer, Serialize};
use std::fmt::Write as StdFmtWrite;

#[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub struct AsciiChar(u8);

impl<'de> Deserialize<'de> for AsciiChar {
  fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
  where
    D: Deserializer<'de>,
  {
    struct AsciiCharVisitor;

    impl Visitor<'_> for AsciiCharVisitor {
      type Value = AsciiChar;

      fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        formatter.write_str("an ASCII byte value (0-127)")
      }

      fn visit_u8<E>(self, value: u8) -> Result<Self::Value, E>
      where
        E: de::Error,
      {
        if value < 128 {
          Ok(AsciiChar(value))
        } else {
          Err(E::custom(format!("value {value} is not ASCII (must be < 128)")))
        }
      }

      fn visit_u64<E>(self, value: u64) -> Result<Self::Value, E>
      where
        E: de::Error,
      {
        if value < 128 {
          Ok(AsciiChar(value as u8))
        } else {
          Err(E::custom(format!("value {value} is not ASCII (must be < 128)")))
        }
      }
    }

    deserializer.deserialize_u8(AsciiCharVisitor)
  }
}

impl AsciiChar {
  /// Create an `AsciiChar` from a byte value.
  ///
  /// # Panics
  ///
  /// Panics if `value` is not ASCII (>= 128).
  pub const fn new(value: u8) -> Self {
    assert!(value < 128, "AsciiChar::new: value >= 128");
    Self(value)
  }

  /// Create an `AsciiChar` from a pre-validated byte value.
  ///
  /// # Precondition
  ///
  /// The caller must ensure that `value` is ASCII (< 128). Passing non-ASCII
  /// input violates the type invariant and causes undefined behavior when
  /// the resulting `AsciiChar` is used in a `Seq` and converted to `&str`.
  ///
  /// Use [`new`](Self::new) or [`From<u8>`] for untrusted input.
  pub const fn from_byte_unchecked(value: u8) -> Self {
    debug_assert!(value < 128, "AsciiChar::from_byte_unchecked: value >= 128");
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

/// # Panics
///
/// Panics if `item` is not ASCII (>= 128).
impl From<u8> for AsciiChar {
  fn from(item: u8) -> Self {
    assert!(item < 128, "AsciiChar::from(u8): value {item} >= 128");
    Self(item)
  }
}

/// # Panics
///
/// Panics if `item` is not ASCII (>= 128).
impl From<u16> for AsciiChar {
  fn from(item: u16) -> Self {
    assert!(item < 128, "AsciiChar::from(u16): value {item} >= 128");
    Self(item as u8)
  }
}

/// # Panics
///
/// Panics if `item` is not ASCII (>= 128).
impl From<u32> for AsciiChar {
  fn from(item: u32) -> Self {
    assert!(item < 128, "AsciiChar::from(u32): value {item} >= 128");
    Self(item as u8)
  }
}

/// # Panics
///
/// Panics if `item` is not ASCII (>= 128).
impl From<u64> for AsciiChar {
  fn from(item: u64) -> Self {
    assert!(item < 128, "AsciiChar::from(u64): value {item} >= 128");
    Self(item as u8)
  }
}

/// # Panics
///
/// Panics if `item` is not ASCII (>= 128).
impl From<usize> for AsciiChar {
  fn from(item: usize) -> Self {
    assert!(item < 128, "AsciiChar::from(usize): value {item} >= 128");
    Self(item as u8)
  }
}

/// # Panics
///
/// Panics if `item` is not ASCII (code point >= 128).
impl From<char> for AsciiChar {
  fn from(item: char) -> Self {
    assert!(item.is_ascii(), "AsciiChar::from(char): '{item}' is not ASCII");
    Self(item as u8)
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
