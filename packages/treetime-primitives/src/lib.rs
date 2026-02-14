mod __tests__;
pub mod bitset128;
pub mod seq;
pub mod seq_char;

pub use bitset128::{BitSet128, Bitset128Status};
pub use seq::Seq;
pub use seq_char::AsciiChar;

/// Minimal trait for alphabet-like types used by I/O operations.
/// Allows treetime-io to validate characters without depending on the full Alphabet type.
pub trait AlphabetLike {
  fn contains(&self, c: AsciiChar) -> bool;
  fn chars(&self) -> impl Iterator<Item = AsciiChar>;
}

impl<A: AlphabetLike> AlphabetLike for &A {
  fn contains(&self, c: AsciiChar) -> bool {
    (*self).contains(c)
  }

  fn chars(&self) -> impl Iterator<Item = AsciiChar> {
    (*self).chars()
  }
}

pub type StateSet = BitSet128;
pub type StateSetStatus = Bitset128Status;

#[macro_export]
macro_rules! stateset {
  ($($args:tt)*) => {
    $crate::bitset128!($($args)*)
  };
}
