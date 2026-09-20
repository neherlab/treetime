#[cfg(test)]
mod tests {
  #![allow(clippy::redundant_clone)]
  use crate::interval::range_complement::range_complement;
  use crate::interval::range_intersection::range_intersection;
  use crate::interval::range_union::range_union;
  use rstest::rstest;
  use std::slice::from_ref;

  #[rstest]
  fn test_range_properties_union_distribution_over_intersection() {
    let set_a = vec![(1, 5), (10, 15)];
    let set_b = vec![(3, 8), (12, 18)];
    let set_c = vec![(6, 10), (14, 20)];

    let actual1 = range_union(&[set_a.clone(), range_intersection(&[set_b.clone(), set_c.clone()])]);
    let actual2 = range_intersection(&[
      range_union(&[set_a.clone(), set_b.clone()]),
      range_union(&[set_a.clone(), set_c.clone()]),
    ]);

    assert_eq!(actual1, actual2);
  }

  #[rstest]
  fn test_range_properties_intersection_absorbs_union() {
    let set_a = vec![(1, 5), (10, 15)];
    let set_b = vec![(3, 8), (12, 18)];

    let actual = range_intersection(&[set_a.clone(), range_union(&[set_a.clone(), set_b.clone()])]);

    assert_eq!(actual, set_a);
  }

  #[rstest]
  fn test_range_properties_union_absorbs_intersection() {
    let set_a = vec![(1, 5), (10, 15)];
    let set_b = vec![(3, 8), (12, 18)];

    let actual = range_union(&[set_a.clone(), range_intersection(&[set_a.clone(), set_b.clone()])]);

    assert_eq!(actual, set_a);
  }

  #[rstest]
  fn test_range_union_with_complement() {
    let set_a = vec![(1, 5), (10, 15)];
    let universe = vec![(0, 20)];

    let actual = range_union(&[set_a.clone(), range_complement(&universe, from_ref(&set_a))]);

    assert_eq!(actual, universe);
  }

  #[rstest]
  fn test_range_union_with_complement_multiple_sets() {
    let set_a = vec![(1, 5), (10, 15)];
    let set_b = vec![(3, 8), (12, 18)];
    let universe = vec![(0, 20)];

    let actual = range_union(&[
      set_a.clone(),
      set_b.clone(),
      range_complement(&universe, &[set_a.clone(), set_b.clone()]),
    ]);

    assert_eq!(actual, universe);
  }

  #[rstest]
  fn test_range_union_complement_empty_set() {
    let universe = vec![(0, 20)];
    let empty_set: Vec<(usize, usize)> = vec![];

    let actual = range_union(&[empty_set.clone(), range_complement(&universe, &[empty_set])]);

    assert_eq!(actual, universe);
  }
}
