use indexmap::IndexSet;
use std::collections::{BTreeSet, HashSet};
use std::hash::Hash;

pub trait SetOperations<T: Eq> {
  #[must_use]
  fn union(&self, other: &Self) -> Self;
  #[must_use]
  fn intersection(&self, other: &Self) -> Self;
}

impl<T: Eq + Hash + Clone> SetOperations<T> for HashSet<T> {
  fn union(&self, other: &Self) -> Self {
    self.union(other).cloned().collect()
  }

  fn intersection(&self, other: &Self) -> Self {
    self.intersection(other).cloned().collect()
  }
}

impl<T: Ord + Clone + Hash> SetOperations<T> for BTreeSet<T> {
  fn union(&self, other: &Self) -> Self {
    self.union(other).cloned().collect()
  }

  fn intersection(&self, other: &Self) -> Self {
    self.intersection(other).cloned().collect()
  }
}

impl<T: Ord + Clone + Hash> SetOperations<T> for IndexSet<T> {
  fn union(&self, other: &Self) -> Self {
    self.union(other).cloned().collect()
  }

  fn intersection(&self, other: &Self) -> Self {
    self.intersection(other).cloned().collect()
  }
}

// Calculate union of multiple sets
pub fn sets_union<T, S, I>(sets: I) -> S
where
  T: Eq + Hash + Clone,
  S: SetOperations<T> + Default,
  I: Iterator<Item = S>,
{
  let mut result = S::default();
  for set in sets {
    result = result.union(&set);
  }
  result
}

// Calculate intersection of multiple sets
pub fn sets_intersection<T, S, I>(sets: I) -> S
where
  T: Eq + Hash + Clone,
  S: SetOperations<T> + Default,
  I: Iterator<Item = S>,
{
  let mut iter = sets;
  let mut result = iter.next().unwrap_or_default();
  for set in iter {
    result = result.intersection(&set);
  }
  result
}

#[cfg(test)]
mod tests {
  #![allow(clippy::redundant_clone)]
  use super::*;
  use maplit::btreeset;

  #[test]
  fn test_sets_union_with_empty() {
    let intersection_set = sets_union(vec![btreeset![1, 2, 5], btreeset![], btreeset![8, 9]].into_iter());
    assert_eq!(intersection_set, btreeset![1, 2, 5, 8, 9]);
  }

  #[test]
  fn test_sets_union_disjoint() {
    let union_set = sets_union(vec![btreeset![1, 2, 3], btreeset![4, 5, 6], btreeset![7, 8, 9]].into_iter());
    assert_eq!(union_set, btreeset![1, 2, 3, 4, 5, 6, 7, 8, 9]);
  }

  #[test]
  fn test_sets_union_overlapping() {
    let union_set = sets_union(vec![btreeset![1, 2, 3], btreeset![3, 4, 5], btreeset![5, 6, 7]].into_iter());
    assert_eq!(union_set, btreeset![1, 2, 3, 4, 5, 6, 7]);
  }

  #[test]
  fn test_sets_union_adjacent() {
    let union_set = sets_union(vec![btreeset![1, 2, 3], btreeset![4, 5, 6], btreeset![6, 7, 8]].into_iter());
    assert_eq!(union_set, btreeset![1, 2, 3, 4, 5, 6, 7, 8]);
  }

  #[test]
  fn test_sets_union_commutativity() {
    let set1 = btreeset![1, 2, 3];
    let set2 = btreeset![2, 3, 4];
    let union_set_1_2 = sets_union(vec![set1.clone(), set2.clone()].into_iter());
    let union_set_2_1 = sets_union(vec![set2, set1].into_iter());
    assert_eq!(union_set_1_2, union_set_2_1);
  }

  #[test]
  fn test_sets_union_associativity() {
    let set1 = btreeset![1, 2, 3];
    let set2 = btreeset![2, 3, 4];
    let set3 = btreeset![3, 4, 5];
    let union_set_1_2_3 = sets_union(vec![set1.clone(), set2.clone(), set3.clone()].into_iter());
    let union_set_2_3_1 = sets_union(vec![set2.clone(), set3.clone(), set1.clone()].into_iter());
    let union_set_3_1_2 = sets_union(vec![set3, set1, set2].into_iter());
    assert_eq!(union_set_1_2_3, union_set_2_3_1);
    assert_eq!(union_set_1_2_3, union_set_3_1_2);
  }

  #[test]
  fn test_sets_union_idempotence() {
    let set1 = btreeset![1, 2, 3];
    let union_set_twice = sets_union(vec![set1.clone(), set1.clone()].into_iter());
    assert_eq!(union_set_twice, set1);
  }

  #[test]
  fn test_sets_intersection_with_empty() {
    let intersection_set = sets_intersection(vec![btreeset![1, 2, 5], btreeset![], btreeset![8, 9]].into_iter());
    assert_eq!(intersection_set, btreeset![]);
  }

  #[test]
  fn test_sets_intersection_disjoint() {
    let intersection_set = sets_intersection(vec![btreeset![1, 2], btreeset![4, 5, 6], btreeset![8, 9]].into_iter());
    assert_eq!(intersection_set, btreeset![]);
  }

  #[test]
  fn test_sets_intersection_overlapping() {
    let intersection_set =
      sets_intersection(vec![btreeset![1, 2, 3], btreeset![2, 3, 4, 5], btreeset![2, 3, 6]].into_iter());
    assert_eq!(intersection_set, btreeset![2, 3]);
  }

  #[test]
  fn test_sets_intersection_adjacent() {
    let intersection_set =
      sets_intersection(vec![btreeset![1, 2, 3], btreeset![4, 5, 6], btreeset![7, 8, 9]].into_iter());
    assert_eq!(intersection_set, btreeset![]);
  }

  #[test]
  fn test_sets_intersection_commutativity() {
    let set1 = btreeset![1, 2, 3];
    let set2 = btreeset![2, 3, 4];
    let intersection_set_1_2 = sets_intersection(vec![set1.clone(), set2.clone()].into_iter());
    let intersection_set_2_1 = sets_intersection(vec![set2, set1].into_iter());
    assert_eq!(intersection_set_1_2, intersection_set_2_1);
  }

  #[test]
  fn test_sets_intersection_associativity() {
    let set1 = btreeset![1, 2, 3];
    let set2 = btreeset![2, 3, 4];
    let set3 = btreeset![3, 4, 5];
    let intersection_set_1_2_3 = sets_intersection(vec![set1.clone(), set2.clone(), set3.clone()].into_iter());
    let intersection_set_2_3_1 = sets_intersection(vec![set2.clone(), set3.clone(), set1.clone()].into_iter());
    let intersection_set_3_1_2 = sets_intersection(vec![set3, set1, set2].into_iter());
    assert_eq!(intersection_set_1_2_3, intersection_set_2_3_1);
    assert_eq!(intersection_set_1_2_3, intersection_set_3_1_2);
  }

  #[test]
  fn test_sets_intersection_idempotence() {
    let set1 = btreeset![1, 2, 3];
    let intersection_set_twice = sets_intersection(vec![set1.clone(), set1.clone()].into_iter());
    assert_eq!(intersection_set_twice, set1);
  }

  // Union distributes over intersection: A ∪ (B ∩ C) = (A ∪ B) ∩ (A ∪ C)
  #[test]
  fn test_sets_union_distributes_over_intersection() {
    let set_a = btreeset![1, 2, 3];
    let set_b = btreeset![2, 3, 4];
    let set_c = btreeset![3, 4, 5];

    // A ∪ (B ∩ C)
    let union_inter = sets_union(
      vec![
        set_a.clone(),
        sets_intersection(vec![set_b.clone(), set_c.clone()].into_iter()),
      ]
      .into_iter(),
    );

    // (A ∪ B) ∩ (A ∪ C)
    let inter_union = sets_intersection(
      vec![
        sets_union(vec![set_a.clone(), set_b].into_iter()),
        sets_union(vec![set_a, set_c].into_iter()),
      ]
      .into_iter(),
    );

    assert_eq!(union_inter, inter_union);
  }

  #[test]
  fn test_sets_intersection_absorbs_union() {
    // Intersection absorbs union: A ∩ (A ∪ B) = A
    let set_a = btreeset![1, 2, 3];
    let set_b = btreeset![2, 3, 4];

    // A ∩ (A ∪ B)
    let inter_union = sets_intersection(
      vec![
        set_a.clone(),
        sets_union(vec![set_a.clone(), set_b.clone()].into_iter()),
      ]
      .into_iter(),
    );

    assert_eq!(inter_union, set_a);
  }

  #[test]
  fn test_sets_union_absorbs_intersection() {
    // Union absorbs intersection: A ∪ (A ∩ B) = A
    let set_a = btreeset![1, 2, 3];
    let set_b = btreeset![2, 3, 4];

    // A ∪ (A ∩ B)
    let union_inter = sets_union(
      vec![
        set_a.clone(),
        sets_intersection(vec![set_a.clone(), set_b.clone()].into_iter()),
      ]
      .into_iter(),
    );

    assert_eq!(union_inter, set_a);
  }
}
