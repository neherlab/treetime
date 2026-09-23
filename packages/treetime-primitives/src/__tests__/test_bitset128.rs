#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::{AsciiChar, BitSet128, BitSet128Status, bitset128};
  use itertools::{Itertools as _, izip};
  use pretty_assertions::{assert_eq, assert_ne};
  use rstest::rstest;
  use std::hash::{DefaultHasher, Hash, Hasher};
  use std::iter::repeat_n;

  #[test]
  fn test_bitset128_new_is_empty() {
    let actual = BitSet128::new();
    let expected = bitset128! {};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_iter() {
    let actual = BitSet128::from_iter(['a', 'b', 'c']);
    let expected = bitset128! {'a', 'b', 'c'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_slice() {
    let actual = BitSet128::from_iter(&['x', 'y', 'z']);
    let expected = bitset128! {'x', 'y', 'z'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_is_empty() {
    let a = BitSet128::new();
    assert!(a.is_empty());
    let a = bitset128! {'a'};
    assert!(!a.is_empty());
  }

  #[test]
  fn test_bitset128_len() {
    let actual = bitset128! {'a', 'b', 'c'};
    let expected_len = 3;
    assert_eq!(actual.len(), expected_len);

    let actual = BitSet128::new();
    let expected_len = 0;
    assert_eq!(actual.len(), expected_len);
  }

  #[test]
  fn test_bitset128_clear() {
    let mut actual = bitset128! {'a', 'b', 'c'};
    actual.clear();
    let expected = BitSet128::new();
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_insert() {
    let mut actual = bitset128! {};
    actual.insert('a');
    actual.insert('a');
    actual.insert('b');
    let expected = bitset128! {'a', 'b'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_remove() {
    let mut actual = bitset128! {'a', 'b', 'c'};
    actual.remove('b');
    actual.remove('b');
    let expected = bitset128! {'a', 'c'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_union() {
    let a = bitset128! {'a', 'b', 'y'};
    let b = bitset128! {'b', 'z', 'x'};
    let actual = a.union(&b);
    let expected = bitset128! {'a', 'b', 'x', 'y', 'z'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_union_with_empty() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {};
    let actual = a.union(&b);
    let expected = bitset128! {'a', 'b', 'c'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_intersection() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {'b', 'c', 'x'};
    let actual = a.intersection(&b);
    let expected = bitset128! {'b', 'c'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_intersection_with_empty() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {};
    let actual = a.intersection(&b);
    let expected = bitset128! {};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_union() {
    let a = bitset128! {'a', 'b', 'x', 'z'};
    let b = bitset128! {'y', 'x', 'a', 'z'};
    let c = bitset128! {'p', 'q', 'y', 'a'};
    let actual = BitSet128::from_union([a, b, c]);
    let expected = bitset128! {'a', 'b', 'p', 'q', 'x', 'y', 'z'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_union_with_empty() {
    let a = bitset128! {'a', 'b', 'x'};
    let b = bitset128! {};
    let c = bitset128! {'p', 'q', 'y'};
    let actual = BitSet128::from_union([a, b, c]);
    let expected = bitset128! {'a', 'b', 'p', 'q', 'x', 'y'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_union_of_all_empty() {
    let a = bitset128! {};
    let b = bitset128! {};
    let c = bitset128! {};
    let actual = BitSet128::from_union([a, b, c]);
    let expected = bitset128! {};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_intersection() {
    let a = bitset128! {'a', 'b', 'x', 'z'};
    let b = bitset128! {'y', 'x', 'a', 'z'};
    let c = bitset128! {'x', 'q', 'y', 'a'};
    let actual = BitSet128::from_intersection([a, b, c]);
    let expected = bitset128! {'a', 'x'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_intersection_with_empty() {
    let a = bitset128! {'a', 'b', 'x'};
    let b = bitset128! {};
    let c = bitset128! {'p', 'q', 'y'};
    let actual = BitSet128::from_intersection([a, b, c]);
    let expected = bitset128! {};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_intersection_of_all_empty() {
    let a = bitset128! {};
    let b = bitset128! {};
    let c = bitset128! {};
    let actual = BitSet128::from_intersection([a, b, c]);
    let expected = bitset128! {};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_plurality_breaks_fitch_union_tie() {
    let sets = [bitset128! {'C'}, bitset128! {'C'}, bitset128! {'A'}];
    let actual = BitSet128::from_plurality(&sets);
    let expected = bitset128! {'C'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_plurality_large_majority() {
    let sets = [
      bitset128! {'C'},
      bitset128! {'C'},
      bitset128! {'C'},
      bitset128! {'C'},
      bitset128! {'A'},
    ];
    let actual = BitSet128::from_plurality(&sets);
    let expected = bitset128! {'C'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_plurality_retains_genuine_tie() {
    let sets = [bitset128! {'C'}, bitset128! {'C'}, bitset128! {'A'}, bitset128! {'A'}];
    let actual = BitSet128::from_plurality(&sets);
    let expected = bitset128! {'A', 'C'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_plurality_overlapping_sets_tie_above_half() {
    let sets = [bitset128! {'A', 'C'}; 5];
    let actual = BitSet128::from_plurality(&sets);
    let expected = bitset128! {'A', 'C'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_plurality_partial_overlap() {
    let sets = [bitset128! {'A', 'C'}, bitset128! {'A', 'C'}, bitset128! {'A', 'G'}];
    let actual = BitSet128::from_plurality(&sets);
    let expected = bitset128! {'A'};
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case(bitset128! {'A', 'C'}, bitset128! {'G', 'T'}, bitset128! {'A', 'C', 'G', 'T'})]
  #[case(bitset128! {'A', 'C'}, bitset128! {'C', 'G'}, bitset128! {'C'})]
  #[case(bitset128! {'A'}, bitset128! {'A'}, bitset128! {'A'})]
  #[case(bitset128! {'A'}, bitset128! {'C'}, bitset128! {'A', 'C'})]
  fn test_bitset128_from_plurality_two_sets_matches_fitch(
    #[case] a: BitSet128,
    #[case] b: BitSet128,
    #[case] expected: BitSet128,
  ) {
    let sets = [a, b];
    let intersection = BitSet128::from_intersection(sets);
    let fitch = if intersection.is_empty() {
      BitSet128::from_union(sets)
    } else {
      intersection
    };
    assert_eq!(BitSet128::from_plurality(&sets), expected);
    assert_eq!(BitSet128::from_plurality(&sets), fitch);
  }

  #[test]
  fn test_bitset128_from_plurality_single_set() {
    let sets = [bitset128! {'A', 'G'}];
    let actual = BitSet128::from_plurality(&sets);
    let expected = bitset128! {'A', 'G'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_plurality_empty_input() {
    let actual = BitSet128::from_plurality(&[]);
    assert_eq!(actual, bitset128! {});
  }

  #[test]
  fn test_bitset128_from_plurality_all_empty_sets() {
    let sets = [bitset128! {}; 3];
    let actual = BitSet128::from_plurality(&sets);
    assert_eq!(actual, bitset128! {});
  }

  #[rstest]
  #[case(3)]
  #[case(4)]
  fn test_bitset128_from_plurality_matches_brute_force_minimum(#[case] n_children: usize) {
    const STATES: [char; 4] = ['A', 'C', 'G', 'T'];

    let subsets: Vec<BitSet128> = (1_u8..16)
      .map(|mask| {
        let mut set = BitSet128::new();
        for (i, state) in STATES.iter().enumerate() {
          if mask & (1 << i) != 0 {
            set.insert(*state);
          }
        }
        set
      })
      .collect();

    let cases = repeat_n(subsets.iter().copied(), n_children)
      .multi_cartesian_product()
      .collect_vec();

    let disagreements = cases
      .iter()
      .filter_map(|sets| {
        let costs = STATES
          .iter()
          .map(|state| sets.iter().filter(|set| !set.contains(*state)).count())
          .collect_vec();
        let min_cost = *costs.iter().min().unwrap();
        let expected: BitSet128 = izip!(&STATES, &costs)
          .filter(|(_, cost)| **cost == min_cost)
          .map(|(state, _)| *state)
          .collect();
        let actual = BitSet128::from_plurality(sets);
        (actual != expected).then(|| format!("sets {sets:?}: expected {expected:?}, actual {actual:?}"))
      })
      .collect_vec();

    assert_eq!(
      cases.len(),
      subsets.len().pow(n_children as u32),
      "enumerated every case"
    );
    assert_eq!(Vec::<String>::new(), disagreements);
  }

  #[test]
  fn test_bitset128_difference() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {'b', 'c', 'x'};
    let actual = a.difference(&b);
    let expected = bitset128! {'a'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_is_disjoint() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {'x', 'y', 'z'};
    assert!(a.is_disjoint(&b));
  }

  #[test]
  fn test_bitset128_is_subset() {
    let a = bitset128! {'a', 'b'};
    let b = bitset128! {'a', 'b', 'c'};
    assert!(a.is_subset(&b));
  }

  #[test]
  fn test_bitset128_is_superset() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {'a', 'b'};
    assert!(a.is_superset(&b));
  }

  #[test]
  fn test_bitset128_display() {
    let a = bitset128! {'a', 'b', 'c'};
    let actual = a.to_string();
    let expected = "{a, b, c}";
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_debug() {
    let a = bitset128! {'a', 'b', 'c'};
    let actual = format!("{a:?}");
    let expected = "{a, b, c}";
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_eq() {
    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {'a', 'b', 'c'};
    let c = bitset128! {'x', 'y', 'z'};
    assert_eq!(a, b);
    assert_ne!(a, c);
  }

  #[test]
  fn test_bitset128_hash() {
    fn calculate_hash<T: Hash>(t: &T) -> u64 {
      let mut hasher = DefaultHasher::new();
      t.hash(&mut hasher);
      hasher.finish()
    }

    let a = bitset128! {'a', 'b', 'c'};
    let b = bitset128! {'a', 'b', 'c'};
    let c = bitset128! {'x', 'y', 'z'};
    assert_eq!(calculate_hash(&a), calculate_hash(&b));
    assert_ne!(calculate_hash(&a), calculate_hash(&c));
  }

  #[test]
  fn test_bitset128_get_empty() {
    let set = BitSet128::new();
    assert!(matches!(set.get(), BitSet128Status::Empty));
  }

  #[test]
  fn test_bitset128_get_unambiguous() {
    let set = BitSet128::from_char('A');
    assert!(matches!(set.get(), BitSet128Status::Unambiguous(c) if c == AsciiChar::from_byte_unchecked(b'A')));
  }

  #[test]
  fn test_bitset128_get_ambiguous() {
    let set = BitSet128::from_iter(vec!['A', 'C']);
    assert!(matches!(set.get(), BitSet128Status::Ambiguous(_)));
    if let BitSet128Status::Ambiguous(actual) = set.get() {
      let expected = bitset128! {'A', 'C'};
      assert_eq!(actual, expected);
    }
  }

  #[test]
  fn test_bitset128_get_one() {
    let set = BitSet128::from_iter(['T', 'A']);
    let actual = set.get_one();
    let expected = AsciiChar::from_byte_unchecked(b'A');
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_get_one_exactly() {
    let set = BitSet128::from_iter(['T']);
    let actual = set.get_one();
    let expected = AsciiChar::from_byte_unchecked(b'T');
    assert_eq!(actual, expected);
  }
}
