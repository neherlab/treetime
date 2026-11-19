#[cfg(test)]
mod tests {
  #![expect(
    clippy::redundant_clone,
    reason = "tests clone inputs to compare them after the call"
  )]
  use crate::interval::range_intersection::range_intersection;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::no_sets(                vec![],                                                                              vec![])]
  #[case::empty_set_annihilates(  vec![vec![(1, 5), (8, 10)], vec![]],                                                 vec![])]
  #[case::single_set(             vec![vec![(1, 5), (8, 10)]],                                                         vec![(1, 5), (8, 10)])]
  #[case::disjoint_sets(          vec![vec![(1, 5)], vec![(8, 10)], vec![(15, 20)]],                                   vec![])]
  #[case::contained_range(        vec![vec![(1, 5)], vec![(2, 4)]],                                                    vec![(2, 4)])]
  #[case::chained_overlap(        vec![vec![(1, 5), (3, 6)], vec![(3, 6), (5, 8)], vec![(4, 7), (6, 9)]],              vec![(4, 6)])]
  #[case::nested_overlap(         vec![vec![(1, 10), (15, 20)], vec![(5, 8), (12, 18)], vec![(2, 6), (16, 22)]],       vec![(5, 6), (16, 18)])]
  #[case::partial_overlap(        vec![vec![(1, 5), (10, 15)], vec![(3, 8), (12, 18)], vec![(6, 10), (14, 20)]],       vec![(14, 15)])]
  #[case::identical_sets(         vec![vec![(1, 5), (8, 10)], vec![(1, 5), (8, 10)], vec![(1, 5), (8, 10)]],           vec![(1, 5), (8, 10)])]
  #[case::shared_start(           vec![vec![(1, 5), (8, 10)], vec![(1, 6), (8, 12)], vec![(1, 7), (8, 14)]],           vec![(1, 5), (8, 10)])]
  #[case::different_range_counts( vec![vec![(1, 5), (8, 10), (15, 20)], vec![(1, 5), (8, 10)]],                       vec![(1, 5), (8, 10)])]
  #[trace]
  fn test_range_intersection_keeps_common_ranges(
    #[case] input: Vec<Vec<(usize, usize)>>,
    #[case] expected: Vec<(usize, usize)>,
  ) {
    let actual = range_intersection(&input);
    assert_eq!(expected, actual);
  }

  #[rstest]
  fn test_range_intersection_commutative() {
    let set1 = vec![(1, 5), (10, 15)];
    let set2 = vec![(3, 8), (12, 18)];
    let actual1 = range_intersection(&[set1.clone(), set2.clone()]);
    let actual2 = range_intersection(&[set2.clone(), set1.clone()]);
    assert_eq!(actual1, actual2);
  }

  #[rstest]
  fn test_range_intersection_associative() {
    let set1 = vec![(1, 5), (10, 15)];
    let set2 = vec![(3, 8), (12, 18)];
    let set3 = vec![(6, 10), (14, 20)];
    let all_at_once = range_intersection(&[set1.clone(), set2.clone(), set3.clone()]);
    let left_first = range_intersection(&[range_intersection(&[set1.clone(), set2.clone()]), set3.clone()]);
    let right_first = range_intersection(&[set1.clone(), range_intersection(&[set2.clone(), set3.clone()])]);
    assert_eq!(all_at_once, left_first);
    assert_eq!(all_at_once, right_first);
  }

  #[rstest]
  fn test_range_intersection_repeated_operand_idempotent() {
    let set1 = vec![(1, 5), (10, 15)];
    let set2 = vec![(3, 8), (12, 18)];
    let actual = range_intersection(&[set1.clone(), set2.clone(), set1.clone()]);
    assert_eq!(range_intersection(&[set1.clone(), set2.clone()]), actual);
  }
}
