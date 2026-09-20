#[cfg(test)]
mod tests {
  #![allow(clippy::redundant_clone)]
  use crate::interval::range_intersection::*;
  use rstest::rstest;

  #[rstest]
  fn test_range_intersection_empty_input() {
    let actual = range_intersection(&[]);
    assert_eq!(actual, vec![]);
  }

  #[rstest]
  fn test_range_intersection_with_empty() {
    let actual = range_intersection(&[vec![(1, 5)], vec![], vec![(15, 20)]]);
    assert_eq!(actual, vec![]);
  }

  #[rstest]
  fn test_range_intersection_single_set() {
    let actual = range_intersection(&[vec![(1, 5), (8, 10)]]);
    assert_eq!(vec![(1, 5), (8, 10)], actual);
  }

  #[rstest]
  fn test_range_intersection_multiple_sets_no_overlap() {
    let actual = range_intersection(&[vec![(1, 5)], vec![(8, 10)], vec![(15, 20)]]);
    assert_eq!(actual, vec![]);
  }

  #[rstest]
  fn test_range_intersection_multiple_sets_with_overlap() {
    let actual = range_intersection(&[vec![(1, 5), (3, 6)], vec![(3, 6), (5, 8)], vec![(4, 7), (6, 9)]]);
    assert_eq!(vec![(4, 6)], actual);
  }

  #[rstest]
  fn test_range_intersection_multiple_sets_with_nested_overlap() {
    let actual = range_intersection(&[vec![(1, 10), (15, 20)], vec![(5, 8), (12, 18)], vec![(2, 6), (16, 22)]]);
    assert_eq!(vec![(5, 6), (16, 18)], actual);
  }

  #[rstest]
  fn test_range_intersection_multiple_sets_with_same_ranges() {
    let actual = range_intersection(&[vec![(1, 5), (8, 10)], vec![(1, 5), (8, 10)], vec![(1, 5), (8, 10)]]);
    assert_eq!(vec![(1, 5), (8, 10)], actual);
  }

  #[rstest]
  fn test_range_intersection_single_interval_contained_within_other() {
    let actual = range_intersection(&[vec![(1, 5)], vec![(2, 4)]]);
    assert_eq!(vec![(2, 4)], actual);
  }

  #[rstest]
  fn test_range_intersection_single_interval_not_contained_within_other() {
    let actual = range_intersection(&[vec![(1, 5)], vec![(6, 8)]]);
    assert_eq!(actual, vec![]);
  }

  #[rstest]
  fn test_range_intersection_multiple_sets_with_partial_overlap() {
    let actual = range_intersection(&[vec![(1, 5), (10, 15)], vec![(3, 8), (12, 18)], vec![(6, 10), (14, 20)]]);
    assert_eq!(vec![(14, 15)], actual);
  }

  #[rstest]
  fn test_range_intersection_sets_with_same_start_or_end() {
    let actual = range_intersection(&[vec![(1, 5), (8, 10)], vec![(1, 6), (8, 12)], vec![(1, 7), (8, 14)]]);
    assert_eq!(vec![(1, 5), (8, 10)], actual);
  }

  #[rstest]
  fn test_range_intersection_sets_with_same_range_different_count() {
    let actual = range_intersection(&[vec![(1, 5), (8, 10), (15, 20)], vec![(1, 5), (8, 10)]]);
    assert_eq!(vec![(1, 5), (8, 10)], actual);
  }

  #[rstest]
  fn test_range_intersection_sets_with_empty_intervals() {
    let actual = range_intersection(&[vec![(1, 5), (8, 10)], vec![], vec![(15, 20)]]);
    assert_eq!(actual, vec![]);
  }

  #[rstest]
  fn test_range_intersection_commutativity() {
    let set1 = vec![(1, 5), (10, 15)];
    let set2 = vec![(3, 8), (12, 18)];
    let actual1 = range_intersection(&[set1.clone(), set2.clone()]);
    let actual2 = range_intersection(&[set2.clone(), set1.clone()]);
    assert_eq!(actual1, actual2);
  }

  #[rstest]
  fn test_range_intersection_associativity() {
    let set1 = vec![(1, 5), (10, 15)];
    let set2 = vec![(3, 8), (12, 18)];
    let set3 = vec![(6, 10), (14, 20)];
    let actual1 = range_intersection(&[set1.clone(), set2.clone(), set3.clone()]);
    let actual2 = range_intersection(&[set1.clone(), set2.clone()]);
    let actual3 = range_intersection(&[set2.clone(), set3.clone()]);
    let actual2_3 = range_intersection(&[actual2.clone(), set3.clone()]);
    assert_eq!(actual1, actual2_3);
  }

  #[rstest]
  fn test_range_intersection_idempotence() {
    let set1 = vec![(1, 5), (10, 15)];
    let actual = range_intersection(&[set1.clone(), set1.clone()]);
    assert_eq!(actual, set1);
  }

  #[rstest]
  fn test_range_intersection_absorption() {
    let set1 = vec![(1, 5), (10, 15)];
    let set2 = vec![(3, 8), (12, 18)];
    let actual = range_intersection(&[set1.clone(), set2.clone(), set1.clone()]);
    assert_eq!(actual, range_intersection(&[set1.clone(), set2.clone()]));
  }

  #[rstest]
  fn test_range_intersection_empty_set() {
    let set1 = vec![(1, 5), (10, 15)];
    let actual = range_intersection(&[set1.clone(), vec![]]);
    assert_eq!(actual, vec![]);
  }
}
