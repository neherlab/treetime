#[cfg(test)]
mod tests {
  #![expect(
    clippy::redundant_clone,
    reason = "tests clone inputs to compare them after the call"
  )]
  use crate::interval::range_union::range_union;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::no_sets(               vec![],                                                                 vec![])]
  #[case::empty_sets(            vec![vec![], vec![]],                                                   vec![])]
  #[case::empty_set_is_identity( vec![vec![(1, 5)], vec![], vec![(15, 20)]],                             vec![(1, 5), (15, 20)])]
  #[case::single_set(            vec![vec![(1, 5), (8, 10)]],                                            vec![(1, 5), (8, 10)])]
  #[case::disjoint_sets(         vec![vec![(1, 5), (8, 10)], vec![(15, 20), (25, 30)]],                  vec![(1, 5), (8, 10), (15, 20), (25, 30)])]
  #[case::partial_overlap(       vec![vec![(1, 5), (8, 10)], vec![(3, 6), (12, 15)]],                    vec![(1, 6), (8, 10), (12, 15)])]
  #[case::chained_overlap(       vec![vec![(1, 5), (3, 6)], vec![(3, 6), (5, 8)], vec![(4, 7), (6, 9)]], vec![(1, 9)])]
  #[case::nested_overlap(        vec![vec![(1, 10), (15, 20)], vec![(5, 8), (12, 18)], vec![(2, 6), (16, 22)]], vec![(1, 10), (12, 22)])]
  #[case::identical_sets(        vec![vec![(1, 5), (8, 10)], vec![(1, 5), (8, 10)], vec![(1, 5), (8, 10)]],     vec![(1, 5), (8, 10)])]
  #[case::subset_absorbed(       vec![vec![(1, 5), (10, 15)], vec![(2, 4)]],                             vec![(1, 5), (10, 15)])]
  #[trace]
  fn test_range_union_merges_ranges(#[case] input: Vec<Vec<(usize, usize)>>, #[case] expected: Vec<(usize, usize)>) {
    let actual = range_union(&input);
    assert_eq!(expected, actual);
  }

  #[rstest]
  fn test_range_union_commutative() {
    let input1 = vec![(1, 5), (10, 15)];
    let input2 = vec![(8, 10), (20, 25)];

    let union1 = range_union(&[input1.clone(), input2.clone()]);
    let union2 = range_union(&[input2, input1]);

    assert_eq!(union1, union2);
  }

  #[rstest]
  fn test_range_union_associative() {
    let input1 = vec![(1, 5), (10, 15)];
    let input2 = vec![(8, 10), (20, 25)];
    let input3 = vec![(18, 20), (30, 35)];

    let union1 = range_union(&[input1.clone(), input2.clone(), input3.clone()]);
    let union2 = range_union(&[input1.clone(), range_union(&[input2, input3])]);

    assert_eq!(union1, union2);
  }
}
