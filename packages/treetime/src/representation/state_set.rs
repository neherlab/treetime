pub use crate::representation::bitset128::BitSet128;
pub use crate::representation::bitset128::Bitset128Status;

pub type StateSetStatus = Bitset128Status;
pub type StateSet = BitSet128;

#[macro_export]
macro_rules! stateset {
  ($($args:tt)*) => {
    $crate::bitset128!($($args)*)
  };
}
