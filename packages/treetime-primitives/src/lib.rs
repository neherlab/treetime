mod bitset128;
mod seq;
mod seq_char;

pub use bitset128::{BitSet128, Bitset128Status};
pub use seq::Seq;
pub use seq_char::AsciiChar;

pub type StateSet = BitSet128;
pub type StateSetStatus = Bitset128Status;

#[macro_export]
macro_rules! stateset {
  ($($args:tt)*) => {
    $crate::bitset128!($($args)*)
  };
}
