#[cfg(test)]
mod __tests__;

use crate::interval::range::{from_interval_set, to_interval_sets};
use gcollections::ops::{Empty, Union};
use intervallum::interval_set::IntervalSet;

pub fn range_union(range_sets: &[Vec<(usize, usize)>]) -> Vec<(usize, usize)> {
  range_union_iter(range_sets.iter()).collect()
}

pub fn range_union_iter<'a>(
  range_sets: impl Iterator<Item = &'a Vec<(usize, usize)>> + 'a,
) -> impl Iterator<Item = (usize, usize)> {
  let iter = to_interval_sets(range_sets);
  let result = interval_sets_union(iter);
  from_interval_set(result)
}

pub fn interval_sets_union(interval_sets: impl Iterator<Item = IntervalSet<usize>>) -> IntervalSet<usize> {
  interval_sets.fold(IntervalSet::empty(), |acc, set| acc.union(&set))
}
