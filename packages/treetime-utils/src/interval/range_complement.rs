use crate::interval::range::{from_interval_set, to_interval_set, to_interval_sets};
use crate::interval::range_union::interval_sets_union;
use gcollections::ops::Difference;
use itertools::Itertools;

pub fn range_complement(universe: &[(usize, usize)], range_sets: &[Vec<(usize, usize)>]) -> Vec<(usize, usize)> {
  range_complement_iter(universe, range_sets.iter()).collect_vec()
}

pub fn range_complement_iter<'a>(
  universe: &[(usize, usize)],
  range_sets: impl Iterator<Item = &'a Vec<(usize, usize)>> + 'a,
) -> impl Iterator<Item = (usize, usize)> {
  let interval_sets = to_interval_sets(range_sets);
  let union = interval_sets_union(interval_sets);
  let universe = to_interval_set(universe);
  let result = universe.difference(&union);
  from_interval_set(result)
}

#[cfg(test)]
mod tests {
  #![allow(clippy::redundant_clone)]
  use super::*;
  use rstest::rstest;

  #[rstest]
  fn test_range_complement_empty() {
    let universe = vec![(0, 10)];
    let interval_sets: Vec<Vec<(usize, usize)>> = vec![];

    let complement_set = range_complement(&universe, &interval_sets);
    assert_eq!(vec![(0, 10)], complement_set);
  }

  #[rstest]
  fn test_range_complement_multiple() {
    let universe = vec![(0, 10)];
    let interval_sets = vec![vec![(2, 4), (7, 9)], vec![(1, 2), (5, 6)]];

    let actual = range_complement(&universe, &interval_sets);
    let expected = vec![(0, 1), (4, 5), (6, 7), (9, 10)];
    assert_eq!(actual, expected);
  }

  #[rstest]
  fn test_range_complement_full_universe() {
    let universe = vec![(0, 10)];
    let interval_sets = vec![vec![(0, 10)]];

    let complement_set = range_complement(&universe, &interval_sets);
    assert_eq!(complement_set, vec![]);
  }

  #[rstest]
  fn test_range_complement_full_universe_overlap() {
    let universe = vec![(0, 10)];
    let interval_sets = vec![vec![(0, 10)], vec![(6, 10)]];

    let complement_set = range_complement(&universe, &interval_sets);
    assert_eq!(complement_set, vec![]);
  }

  #[rstest]
  fn test_range_complement_full_universe_adjacent() {
    let universe = vec![(0, 10)];
    let interval_sets = vec![vec![(0, 3)], vec![(3, 10)]];

    let complement_set = range_complement(&universe, &interval_sets);
    assert_eq!(complement_set, vec![]);
  }

  #[rstest]
  fn test_range_complement_full_universe_same() {
    let universe = vec![(0, 10)];
    let interval_sets = vec![vec![(0, 10)], vec![(0, 10)]];

    let complement_set = range_complement(&universe, &interval_sets);
    assert_eq!(complement_set, vec![]);
  }
}
