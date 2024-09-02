use crate::utils::interval::range::{from_interval_set, to_interval_sets};
use gcollections::ops::{Empty, Union};
use interval::interval_set::IntervalSet;

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

#[cfg(test)]
mod tests {
  #![allow(clippy::redundant_clone)]
  use super::*;
  use rstest::rstest;

  #[rstest]
  fn test_range_union_empty() {
    let actual = range_union(&[]);
    assert_eq!(actual, vec![]);
  }

  #[rstest]
  fn test_range_union_empty_multiple() {
    let actual = range_union(&[vec![], vec![]]);
    assert_eq!(actual, vec![]);
  }

  #[rstest]
  fn test_range_union_with_empty() {
    let actual = range_union(&[
      vec![(1, 5)],
      vec![], // Union with empty is the same as the non-empty set
      vec![(15, 20)],
    ]);
    assert_eq!(actual, vec![(1, 5), (15, 20)]);
  }

  #[rstest]
  fn test_range_union_single_set() {
    let actual = range_union(&[vec![(1, 5), (8, 10)]]);
    assert_eq!(actual, vec![(1, 5), (8, 10)]);
  }

  #[rstest]
  fn test_range_union_multiple_sets_no_overlap() {
    let actual = range_union(&[
      vec![(1, 5)],   // No overlap with other sets
      vec![(8, 10)],  // No overlap with other sets
      vec![(15, 20)], // No overlap with other sets
    ]);
    assert_eq!(actual, vec![(1, 5), (8, 10), (15, 20)]);
  }

  #[rstest]
  fn test_range_union_multiple_sets_with_overlap() {
    let actual = range_union(&[
      vec![(1, 5), (3, 6)], // Overlaps with the second set
      vec![(3, 6), (5, 8)], // Overlaps with the first and third set
      vec![(4, 7), (6, 9)], // Overlaps with the first and second set
    ]);
    assert_eq!(actual, vec![(1, 9)]);
  }

  #[rstest]
  fn test_range_union_multiple_sets_with_nested_overlap() {
    let actual = range_union(&[
      vec![(1, 10), (15, 20)], // No overlaps
      vec![(5, 8), (12, 18)],  // Overlaps with Set 1
      vec![(2, 6), (16, 22)],  // Overlaps with Set 1 and Set 2
    ]);
    assert_eq!(actual, vec![(1, 10), (12, 22)]);
  }

  #[rstest]
  fn test_range_union_multiple_sets_with_same_ranges() {
    let actual = range_union(&[vec![(1, 5), (8, 10)], vec![(1, 5), (8, 10)], vec![(1, 5), (8, 10)]]);
    assert_eq!(actual, vec![(1, 5), (8, 10)]);
  }

  #[rstest]
  fn test_range_union_disjoint_sets() {
    let actual = range_union(&[
      vec![(1, 5), (8, 10)],    //
      vec![(15, 20), (25, 30)], //
    ]);
    assert_eq!(actual, vec![(1, 5), (8, 10), (15, 20), (25, 30)]);
  }

  #[rstest]
  fn test_range_union_overlapping_sets_different_lengths() {
    let actual = range_union(&[
      vec![(1, 5), (8, 10)],  //
      vec![(3, 6), (12, 15)], //
    ]);
    assert_eq!(actual, vec![(1, 6), (8, 10), (12, 15)]);
  }

  #[rstest]
  fn test_range_union_commutativity() {
    let input1 = vec![(1, 5), (10, 15)];
    let input2 = vec![(8, 10), (20, 25)];

    let union1 = range_union(&[input1.clone(), input2.clone()]);
    let union2 = range_union(&[input2, input1]);

    assert_eq!(union1, union2);
  }

  #[rstest]
  fn test_range_union_associativity() {
    let input1 = vec![(1, 5), (10, 15)];
    let input2 = vec![(8, 10), (20, 25)];
    let input3 = vec![(18, 20), (30, 35)];

    let union1 = range_union(&[input1.clone(), input2.clone(), input3.clone()]);
    let union2 = range_union(&[input1.clone(), range_union(&[input2, input3])]);

    assert_eq!(union1, union2);
  }

  #[rstest]
  fn test_range_union_idempotence() {
    let input = vec![(1, 5), (10, 15)];

    let union1 = range_union(&[input.clone(), input.clone()]);
    let union2 = range_union(&[input]);

    assert_eq!(union1, union2);
  }

  #[rstest]
  fn test_range_union_empty_set_identity() {
    let input = vec![(1, 5), (10, 15)];

    let union1 = range_union(&[input.clone(), vec![]]);
    let union2 = range_union(&[input]);

    assert_eq!(union1, union2);
  }

  #[rstest]
  fn test_range_union_absorption() {
    let input1 = vec![(1, 5), (10, 15)];
    let input2 = vec![(2, 4)];
    let union_result = range_union(&[input1.clone(), input2.clone()]);
    // Since input2 is fully encompassed by input1, the union should be equal to input1
    assert_eq!(union_result, input1);
  }
}
