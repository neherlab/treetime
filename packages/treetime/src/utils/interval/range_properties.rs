#[cfg(test)]
mod tests {
  #![allow(clippy::redundant_clone)]
  use crate::utils::interval::range_complement::range_complement;
  use crate::utils::interval::range_intersection::range_intersection;
  use crate::utils::interval::range_union::range_union;
  use rstest::rstest;

  #[rstest]
  fn test_range_properties_union_distribution_over_intersection() {
    // Test for distributivity property:
    // Union distributes over intersection: A ∪ (B ∩ C) = (A ∪ B) ∩ (A ∪ C)

    let set_a = vec![(1, 5), (10, 15)]; // Set A
    let set_b = vec![(3, 8), (12, 18)]; // Set B
    let set_c = vec![(6, 10), (14, 20)]; // Set C

    // A ∪ (B ∩ C)
    let actual1 = range_union(&[set_a.clone(), range_intersection(&[set_b.clone(), set_c.clone()])]);
    // (A ∪ B) ∩ (A ∪ C)
    let actual2 = range_intersection(&[
      range_union(&[set_a.clone(), set_b.clone()]),
      range_union(&[set_a.clone(), set_c.clone()]),
    ]);

    assert_eq!(actual1, actual2);
  }

  #[rstest]
  fn test_range_properties_intersection_absorbs_union() {
    // Test for absorption property:
    // Intersection absorbs union: A ∩ (A ∪ B) = A

    let set_a = vec![(1, 5), (10, 15)]; // Set A
    let set_b = vec![(3, 8), (12, 18)]; // Set B

    // A ∩ (A ∪ B)
    let actual = range_intersection(&[set_a.clone(), range_union(&[set_a.clone(), set_b.clone()])]);

    assert_eq!(actual, set_a);
  }

  #[rstest]
  fn test_range_properties_union_absorbs_intersection() {
    // Test for absorption property in reverse:
    // Union absorbs intersection: A ∪ (A ∩ B) = A

    let set_a = vec![(1, 5), (10, 15)]; // Set A
    let set_b = vec![(3, 8), (12, 18)]; // Set B

    // A ∪ (A ∩ B)
    let actual = range_union(&[set_a.clone(), range_intersection(&[set_a.clone(), set_b.clone()])]);

    assert_eq!(actual, set_a);
  }

  // Test for complementarity property:
  // A ∪ ¬A = U (the universal set), where ¬A represents the complement of A.
  #[rstest]
  fn test_range_union_with_complement() {
    let set_a = vec![(1, 5), (10, 15)]; // Set A
    let universe = vec![(0, 20)]; // Universal set

    // A ∪ ¬A
    let actual = range_union(&[set_a.clone(), range_complement(&universe, &[set_a.clone()])]);

    assert_eq!(actual, universe);
  }

  // Additional test for complementarity property with multiple sets:
  // A ∪ B ∪ ¬(A ∪ B) = U
  #[rstest]
  fn test_range_union_with_complement_multiple_sets() {
    let set_a = vec![(1, 5), (10, 15)]; // Set A
    let set_b = vec![(3, 8), (12, 18)]; // Set B
    let universe = vec![(0, 20)]; // Universal set

    // A ∪ B ∪ ¬(A ∪ B)
    let actual = range_union(&[
      set_a.clone(),
      set_b.clone(),
      range_complement(&universe, &[set_a.clone(), set_b.clone()]),
    ]);

    assert_eq!(actual, universe);
  }

  // Test for complement of empty set:
  // Complement of an empty set within the universe should be the universe itself.
  #[rstest]
  fn test_range_union_complement_empty_set() {
    let universe = vec![(0, 20)];
    let empty_set: Vec<(usize, usize)> = vec![];

    // ∅ ∪ ¬∅
    let actual = range_union(&[empty_set.clone(), range_complement(&universe, &[empty_set])]);

    assert_eq!(actual, universe);
  }
}
