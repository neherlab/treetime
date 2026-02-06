use itertools::{Itertools, chain};

/// Creates a union of two iterators
pub fn iterator_union<I1, I2, T>(iter1: I1, iter2: I2) -> impl Iterator<Item = T>
where
  I1: IntoIterator<Item = T>,
  I2: IntoIterator<Item = T>,
  T: Ord,
{
  chain!(iter1, iter2).sorted().dedup()
}
