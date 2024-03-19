use crate::seq::range::{from_interval_set, to_interval_sets};
use gcollections::ops::{Empty, Intersection};
use interval::interval_set::IntervalSet;

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

#[cfg(test)]
mod tests {
  #![allow(clippy::redundant_clone)]
  use super::*;
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
    let actual = range_intersection(&[
      vec![(1, 5), (8, 10)], // Same range
      vec![(1, 5), (8, 10)], // Same range
      vec![(1, 5), (8, 10)], // Same range
    ]);
    assert_eq!(actual, vec![(1, 5), (8, 10)]);
  }

  #[rstest]
  fn test_range_intersection_single_interval_contained_within_other() {
    let actual = range_intersection(&[vec![(1, 5)], vec![(2, 4)]]);
    assert_eq!(actual, vec![(2, 4)]);
  }

  #[rstest]
  fn test_range_intersection_single_interval_not_contained_within_other() {
    let actual = range_intersection(&[vec![(1, 5)], vec![(6, 8)]]);
    assert_eq!(actual, vec![]);
  }

  #[rstest]
  fn test_range_intersection_multiple_sets_with_partial_overlap() {
    let actual = range_intersection(&[
      vec![(1, 5), (10, 15)],  // No overlap with other sets
      vec![(3, 8), (12, 18)],  // Partial overlap with Set 1
      vec![(6, 10), (14, 20)], // Partial overlap with Set 1 and Set 2
    ]);
    assert_eq!(actual, vec![(14, 15)]);
  }

  #[rstest]
  fn test_range_intersection_sets_with_same_start_or_end() {
    let actual = range_intersection(&[
      vec![(1, 5), (8, 10)], //
      vec![(1, 6), (8, 12)], //
      vec![(1, 7), (8, 14)], //
    ]);
    assert_eq!(actual, vec![(1, 5), (8, 10)]);
  }

  #[rstest]
  fn test_range_intersection_sets_with_same_range_different_count() {
    let actual = range_intersection(&[
      vec![(1, 5), (8, 10), (15, 20)], //
      vec![(1, 5), (8, 10)],           //
    ]);
    assert_eq!(actual, vec![(1, 5), (8, 10)]);
  }

  #[rstest]
  fn test_range_intersection_sets_with_empty_intervals() {
    let actual = range_intersection(&[vec![(1, 5), (8, 10)], vec![], vec![(15, 20)]]);
    assert_eq!(actual, vec![]);
  }

  // Test for commutativity property:
  // Intersection of two sets should yield the same result regardless of the order of sets.
  #[rstest]
  fn test_range_intersection_commutativity() {
    let set1 = vec![(1, 5), (10, 15)]; // Set A
    let set2 = vec![(3, 8), (12, 18)]; // Set B
    let actual1 = range_intersection(&[set1.clone(), set2.clone()]); // A ∩ B
    let actual2 = range_intersection(&[set2.clone(), set1.clone()]); // B ∩ A
    assert_eq!(actual1, actual2);
  }

  #[rstest]
  fn test_range_intersection_associativity() {
    // Test for associativity property:
    // The order of operations should not affect the result.
    let set1 = vec![(1, 5), (10, 15)]; // Set A
    let set2 = vec![(3, 8), (12, 18)]; // Set B
    let set3 = vec![(6, 10), (14, 20)]; // Set C
    let actual1 = range_intersection(&[set1.clone(), set2.clone(), set3.clone()]); // (A ∩ B) ∩ C
    let actual2 = range_intersection(&[set1.clone(), set2.clone()]); // A ∩ B
    let actual3 = range_intersection(&[set2.clone(), set3.clone()]); // B ∩ C
    let actual2_3 = range_intersection(&[actual2.clone(), set3.clone()]); // (A ∩ B) ∩ C
    assert_eq!(actual1, actual2_3);
  }

  #[rstest]
  fn test_range_intersection_idempotence() {
    // Test for idempotence property:
    // Intersecting a set with itself multiple times should not change the result.
    let set1 = vec![(1, 5), (10, 15)]; // Set A
    let actual = range_intersection(&[set1.clone(), set1.clone()]); // A ∩ A
    assert_eq!(actual, set1);
  }

  #[rstest]
  fn test_range_intersection_absorption() {
    // Test for absorption property:
    // If a set is already included in another set, the intersection with the other set should yield the same set.
    let set1 = vec![(1, 5), (10, 15)]; // Set A
    let set2 = vec![(3, 8), (12, 18)]; // Set B
    let actual = range_intersection(&[set1.clone(), set2.clone(), set1.clone()]); // (A ∩ B) ∩ A
    assert_eq!(actual, range_intersection(&[set1.clone(), set2.clone()])); // A ∩ B
  }

  #[rstest]
  fn test_range_intersection_empty_set() {
    // Test for empty set property:
    // Intersecting any set with an empty set should always result in an empty set.
    let set1 = vec![(1, 5), (10, 15)]; // Set A
    let actual = range_intersection(&[set1.clone(), vec![]]); // A ∩ ∅
    assert_eq!(actual, vec![]);
  }
}
