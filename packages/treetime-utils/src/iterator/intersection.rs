use itertools::{EitherOrBoth, Itertools};

/// Creates an intersection of two iterators (items present in both).
///
/// Both inputs are sorted and deduplicated before intersecting.
/// Output is sorted and deduplicated.
pub fn iterator_intersection<I1, I2, T>(iter1: I1, iter2: I2) -> impl Iterator<Item = T>
where
  I1: IntoIterator<Item = T>,
  I2: IntoIterator<Item = T>,
  T: Ord,
{
  iter1
    .into_iter()
    .sorted()
    .dedup()
    .merge_join_by(iter2.into_iter().sorted().dedup(), |a, b| a.cmp(b))
    .filter_map(|pair| match pair {
      EitherOrBoth::Both(item, _) => Some(item),
      _ => None,
    })
}

#[cfg(test)]
mod tests {
  use super::*;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::both_empty(          vec![],        vec![],        vec![]       )]
  #[case::left_empty(          vec![],        vec![1, 2],    vec![]       )]
  #[case::right_empty(         vec![1, 2],    vec![],        vec![]       )]
  #[case::disjoint(            vec![1, 3],    vec![2, 4],    vec![]       )]
  #[case::full_overlap(        vec![1, 2, 3], vec![1, 2, 3], vec![1, 2, 3])]
  #[case::partial_overlap(     vec![1, 2, 3], vec![2, 3, 4], vec![2, 3]   )]
  #[case::single_common(       vec![5],       vec![5],       vec![5]      )]
  #[case::duplicates_in_left(  vec![1, 1, 2], vec![1, 2],    vec![1, 2]   )]
  #[case::duplicates_in_right( vec![1, 2],    vec![2, 2, 3], vec![2]      )]
  #[case::unsorted_input(      vec![3, 1, 2], vec![2, 3, 5], vec![2, 3]   )]
  #[trace]
  fn test_iterator_intersection(#[case] a: Vec<i32>, #[case] b: Vec<i32>, #[case] expected: Vec<i32>) {
    let result: Vec<i32> = iterator_intersection(a, b).collect_vec();
    assert_eq!(expected, result);
  }
}
