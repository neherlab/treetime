#[cfg(test)]
mod tests {
  #![expect(
    clippy::redundant_clone,
    reason = "tests clone inputs to compare them after the call"
  )]
  use crate::interval::range_union::*;
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
    let actual = range_union(&[vec![(1, 5)], vec![], vec![(15, 20)]]);
    assert_eq!(vec![(1, 5), (15, 20)], actual);
  }

  #[rstest]
  fn test_range_union_single_set() {
    let actual = range_union(&[vec![(1, 5), (8, 10)]]);
    assert_eq!(vec![(1, 5), (8, 10)], actual);
  }

  #[rstest]
  fn test_range_union_multiple_sets_no_overlap() {
    let actual = range_union(&[vec![(1, 5)], vec![(8, 10)], vec![(15, 20)]]);
    assert_eq!(vec![(1, 5), (8, 10), (15, 20)], actual);
  }

  #[rstest]
  fn test_range_union_multiple_sets_with_overlap() {
    let actual = range_union(&[vec![(1, 5), (3, 6)], vec![(3, 6), (5, 8)], vec![(4, 7), (6, 9)]]);
    assert_eq!(vec![(1, 9)], actual);
  }

  #[rstest]
  fn test_range_union_multiple_sets_with_nested_overlap() {
    let actual = range_union(&[vec![(1, 10), (15, 20)], vec![(5, 8), (12, 18)], vec![(2, 6), (16, 22)]]);
    assert_eq!(vec![(1, 10), (12, 22)], actual);
  }

  #[rstest]
  fn test_range_union_multiple_sets_with_same_ranges() {
    let actual = range_union(&[vec![(1, 5), (8, 10)], vec![(1, 5), (8, 10)], vec![(1, 5), (8, 10)]]);
    assert_eq!(vec![(1, 5), (8, 10)], actual);
  }

  #[rstest]
  fn test_range_union_disjoint_sets() {
    let actual = range_union(&[vec![(1, 5), (8, 10)], vec![(15, 20), (25, 30)]]);
    assert_eq!(vec![(1, 5), (8, 10), (15, 20), (25, 30)], actual);
  }

  #[rstest]
  fn test_range_union_overlapping_sets_different_lengths() {
    let actual = range_union(&[vec![(1, 5), (8, 10)], vec![(3, 6), (12, 15)]]);
    assert_eq!(vec![(1, 6), (8, 10), (12, 15)], actual);
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
    assert_eq!(union_result, input1);
  }
}
