#[cfg(test)]
mod __tests__;

use crate::interval::range::{from_interval_set, to_interval_sets};
use gcollections::ops::{Empty, Intersection};
use intervallum::interval_set::IntervalSet;

pub fn range_intersection(range_sets: &[Vec<(usize, usize)>]) -> Vec<(usize, usize)> {
  range_intersection_iter(range_sets.iter()).collect()
}

pub fn range_intersection_iter<'a>(
  range_sets: impl Iterator<Item = &'a Vec<(usize, usize)>> + 'a,
) -> impl Iterator<Item = (usize, usize)> {
  let iter = to_interval_sets(range_sets);
  let result = interval_sets_intersection(iter);
  from_interval_set(result)
}

pub fn interval_sets_intersection(mut interval_sets: impl Iterator<Item = IntervalSet<usize>>) -> IntervalSet<usize> {
  let first = interval_sets.next().unwrap_or_else(IntervalSet::empty);
  interval_sets.fold(first, |acc, set| acc.intersection(&set))
}
