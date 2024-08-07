use crate::seq::range::{from_interval_set, to_interval_set};
use gcollections::ops::Difference;

pub fn range_difference(x: &[(usize, usize)], y: &[(usize, usize)]) -> Vec<(usize, usize)> {
  let x = to_interval_set(x);
  let y = to_interval_set(y);
  let result = x.difference(&y);
  from_interval_set(result).collect()
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_range_difference_with_no_overlap() {
    let x = vec![(1, 6), (10, 16)];
    let y = vec![(20, 26)];
    let result = range_difference(&x, &y);
    assert_eq!(result, vec![(1, 6), (10, 16)]);
  }

  #[test]
  fn test_range_difference_with_partial_overlap() {
    let x = vec![(1, 11)];
    let y = vec![(5, 16)];
    let result = range_difference(&x, &y);
    assert_eq!(result, vec![(1, 5)]);
  }

  #[test]
  fn test_range_difference_with_full_overlap() {
    let x = vec![(1, 6)];
    let y = vec![(1, 6)];
    let result = range_difference(&x, &y);
    assert!(result.is_empty());
  }

  #[test]
  fn test_range_difference_with_multiple_overlaps() {
    let x = vec![(1, 6), (6, 11)];
    let y = vec![(2, 8)];
    let result = range_difference(&x, &y);
    assert_eq!(result, vec![(1, 2), (8, 11)]);
  }

  #[test]
  fn test_range_difference_with_sub_intervals() {
    let x = vec![(1, 21)];
    let y = vec![(3, 6), (7, 10), (11, 14), (15, 18)];
    let result = range_difference(&x, &y);
    assert_eq!(result, vec![(1, 3), (6, 7), (10, 11), (14, 15), (18, 21)]);
  }

  #[test]
  fn test_range_difference_with_adjacent_intervals() {
    let x = vec![(1, 6), (6, 11)];
    let y = vec![(2, 7)];
    let result = range_difference(&x, &y);
    assert_eq!(result, vec![(1, 2), (7, 11)]);
  }

  #[test]
  fn test_range_difference_idempotence() {
    let x = vec![(1, 5)];
    let result = range_difference(&x, &x);
    assert!(result.is_empty());
  }

  #[test]
  fn test_range_difference_identity() {
    let x = vec![(1, 5)];
    let empty: Vec<(usize, usize)> = vec![];
    let result = range_difference(&x, &empty);
    assert_eq!(result, x);
  }
}
