use gcollections::ops::{Bounded, Empty, Intersection};
use interval::interval_set::{IntervalSet, ToIntervalSet};

pub fn range_intersection(range_sets: &[Vec<(usize, usize)>]) -> Vec<(usize, usize)> {
  range_intersection_iter(range_sets.iter())
}

#[allow(single_use_lifetimes)]
pub fn range_intersection_iter<'a>(range_sets: impl Iterator<Item = &'a Vec<(usize, usize)>>) -> Vec<(usize, usize)> {
  // Convert input vectors into IntervalSets
  let mut iter = range_sets.map(|a| a.clone().to_interval_set());

  // Perform intersection of all IntervalSets
  let result_set = iter.next().unwrap_or_else(IntervalSet::empty);
  let result_set = iter.fold(result_set, |acc, set| acc.intersection(&set));

  // Convert result set back into vector of tuples
  result_set
    .into_iter()
    .map(|interval| (interval.lower(), interval.upper()))
    .collect::<Vec<_>>()
}

pub fn ranges_contain(ranges: &[(usize, usize)], pos: usize) -> bool {
  ranges.iter().any(|&(start, end)| start <= pos && pos < end)
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_range_intersection_empty_input() {
    let actual = range_intersection(&[]);
    assert_eq!(actual, vec![]);
  }

  #[rstest]
  fn test_range_intersection_with_empty() {
    let actual = range_intersection(&[
      vec![(1, 5)],
      vec![], // Intersection with empty is always empty
      vec![(15, 20)],
    ]);
    assert_eq!(actual, vec![]);
  }

  #[rstest]
  fn test_range_intersection_single_set() {
    let actual = range_intersection(&[vec![(1, 5), (8, 10)]]);
    assert_eq!(actual, vec![(1, 5), (8, 10)]);
  }

  #[rstest]
  fn test_range_intersection_multiple_sets_no_overlap() {
    let actual = range_intersection(&[
      vec![(1, 5)],   // No overlap with other sets
      vec![(8, 10)],  // No overlap with other sets
      vec![(15, 20)], // No overlap with other sets
    ]);
    assert_eq!(actual, vec![]);
  }

  #[rstest]
  fn test_range_intersection_multiple_sets_with_overlap() {
    let actual = range_intersection(&[
      vec![(1, 5), (3, 6)], // Overlaps with the second set
      vec![(3, 6), (5, 8)], // Overlaps with the first and third set
      vec![(4, 7), (6, 9)], // Overlaps with the first and second set
    ]);
    assert_eq!(actual, vec![(4, 6)]);
  }

  #[rstest]
  fn test_range_intersection_multiple_sets_with_nested_overlap() {
    let actual = range_intersection(&[
      vec![(1, 10), (15, 20)], // No overlaps
      vec![(5, 8), (12, 18)],  // Overlaps with Set 1
      vec![(2, 6), (16, 22)],  // Overlaps with Set 1 and Set 2
    ]);
    assert_eq!(actual, vec![(5, 6), (16, 18)]);
  }

  #[rstest]
  fn test_range_intersection_multiple_sets_with_same_ranges() {
    let actual = range_intersection(&[vec![(1, 5), (8, 10)], vec![(1, 5), (8, 10)], vec![(1, 5), (8, 10)]]);
    assert_eq!(actual, vec![(1, 5), (8, 10)]);
  }
}
