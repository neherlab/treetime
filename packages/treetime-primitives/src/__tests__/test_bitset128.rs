#![allow(clippy::op_ref)]

#[cfg(test)]
mod tests {
  use crate::{AsciiChar, BitSet128, BitSet128Status, bitset128};
  use pretty_assertions::{assert_eq, assert_ne};
  use rstest::rstest;
  use std::hash::{DefaultHasher, Hash, Hasher};

  #[test]
  fn test_bitset128_new() {
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
    // The canonical counterexample: Fitch's union fallback yields {A,C} and admits root A at a
    // cost of two changes, while only C is minimal at one.
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
    // Overlapping sets let counts sum above the number of sets, so both states reach count = 5.
    // A simple majority over half the sets therefore cannot imply a unique result.
    let sets = [bitset128! {'A', 'C'}; 5];
    let actual = BitSet128::from_plurality(&sets);
    let expected = bitset128! {'A', 'C'};
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_bitset128_from_plurality_partial_overlap() {
    // count(A) = 3, count(C) = 2, count(G) = 1.
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
    // With one or two sets the plurality coincides with intersect-or-unite, so the existing
    // bifurcating behaviour is preserved exactly.
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

  /// Exhaustive check against a brute-force minimisation of the parsimony cost.
  ///
  /// For a star tree whose children have optimal sets `S_i`, the cost of assigning state `x` to
  /// the parent is `Σ_i [x ∉ S_i]`. This enumerates every combination of non-empty subsets of a
  /// four-letter alphabet, computes the true minimum-cost state set directly from that
  /// definition, and requires `from_plurality` to return exactly it.
  #[rstest]
  #[case(3)]
  #[case(4)]
  fn test_bitset128_from_plurality_matches_brute_force_minimum(#[case] n_children: usize) {
    const STATES: [char; 4] = ['A', 'C', 'G', 'T'];

    // All 15 non-empty subsets of a four-letter alphabet.
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

    let mut n_cases = 0_usize;
    let mut indices = vec![0_usize; n_children];
    loop {
      let sets: Vec<BitSet128> = indices.iter().map(|&i| subsets[i]).collect();

      // Brute force: cost of each state, then every state attaining the minimum.
      let costs: Vec<usize> = STATES
        .iter()
        .map(|state| sets.iter().filter(|set| !set.contains(*state)).count())
        .collect();
      let min_cost = *costs.iter().min().unwrap();
      let mut expected = BitSet128::new();
      for (i, state) in STATES.iter().enumerate() {
        if costs[i] == min_cost {
          expected.insert(*state);
        }
      }

      assert_eq!(
        BitSet128::from_plurality(&sets),
        expected,
        "plurality disagrees with brute-force minimum for sets {sets:?}"
      );
      n_cases += 1;

      // Odometer over subset indices.
      let mut d = 0;
      while d < n_children {
        indices[d] += 1;
        if indices[d] < subsets.len() {
          break;
        }
        indices[d] = 0;
        d += 1;
      }
      if d == n_children {
        break;
      }
    }

    assert_eq!(n_cases, subsets.len().pow(n_children as u32), "enumerated every case");
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

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::disjoint_multi(bitset128!{'a'},            bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case::overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q', 'x'})]
  #[case::identical(bitset128!{'z'},            bitset128!{'z'},            bitset128!{'z'})]
  #[case::partial_overlap(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case::both_empty(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case::multi_disjoint(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e', 'f'},  bitset128!{'a', 'b', 'c', 'd', 'e', 'f'})]
  #[trace]
  fn test_bitset128_add(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a + b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::disjoint_multi(bitset128!{'a'},            bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case::overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q', 'x'})]
  #[case::identical(bitset128!{'z'},            bitset128!{'z'},            bitset128!{'z'})]
  #[case::partial_overlap(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case::both_empty(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case::multi_disjoint(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e', 'f'},  bitset128!{'a', 'b', 'c', 'd', 'e', 'f'})]
  #[trace]
  fn test_bitset128_add_ref_right(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a + &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::disjoint_multi(bitset128!{'a'},            bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case::overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q', 'x'})]
  #[case::identical(bitset128!{'z'},            bitset128!{'z'},            bitset128!{'z'})]
  #[case::partial_overlap(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case::both_empty(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case::multi_disjoint(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e', 'f'},  bitset128!{'a', 'b', 'c', 'd', 'e', 'f'})]
  #[trace]
  fn test_bitset128_add_ref_left(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a + b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::disjoint_multi(bitset128!{'a'},            bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case::overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q', 'x'})]
  #[case::identical(bitset128!{'z'},            bitset128!{'z'},            bitset128!{'z'})]
  #[case::partial_overlap(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case::both_empty(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case::multi_disjoint(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e', 'f'},  bitset128!{'a', 'b', 'c', 'd', 'e', 'f'})]
  #[trace]
  fn test_bitset128_add_ref_both(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a + &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::disjoint_multi(bitset128!{'a'},            bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case::overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q', 'x'})]
  #[case::identical(bitset128!{'z'},            bitset128!{'z'},            bitset128!{'z'})]
  #[case::partial_overlap(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case::both_empty(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case::multi_disjoint(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e', 'f'},  bitset128!{'a', 'b', 'c', 'd', 'e', 'f'})]
  #[trace]
  fn test_bitset128_add_assign(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a += b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::disjoint_multi(bitset128!{'a'},            bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case::overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q', 'x'})]
  #[case::identical(bitset128!{'z'},            bitset128!{'z'},            bitset128!{'z'})]
  #[case::partial_overlap(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case::both_empty(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case::multi_disjoint(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e', 'f'},  bitset128!{'a', 'b', 'c', 'd', 'e', 'f'})]
  #[trace]
  fn test_bitset128_add_assign_ref(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a += &b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case::add_new(bitset128!{'a'},            b'b',                        bitset128!{'a', 'b'})]
  #[case::add_multi(bitset128!{'x', 'y'},       b'z',                        bitset128!{'x', 'y', 'z'})]
  #[case::add_existing(bitset128!{'m', 'n'},       b'm',                        bitset128!{'m', 'n'})]
  #[case::add_to_empty(bitset128!{},               b'a',                        bitset128!{'a'})]
  #[trace]
  fn test_bitset128_add_char_to_set(#[case] a: BitSet128, #[case] b: u8, #[case] expected: BitSet128) {
    let actual = a + b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::add_new(b'a',                        bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::add_multi(b'x',                        bitset128!{'y', 'z'},       bitset128!{'x', 'y', 'z'})]
  #[case::add_existing(b'm',                        bitset128!{'m', 'n'},       bitset128!{'m', 'n'})]
  #[case::add_to_empty(b'a',                        bitset128!{},               bitset128!{'a'})]
  #[trace]
  fn test_bitset128_add_set_to_char(#[case] a: u8, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = b + a;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::both_empty(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'x', 'q'})]
  #[case::overlap_multi(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case::disjoint_multi(bitset128!{'x', 'y'},       bitset128!{'z'},            bitset128!{'x', 'y', 'z'})]
  #[trace]
  fn test_bitset128_or(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a | b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::both_empty(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'x', 'q'})]
  #[case::overlap_multi(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case::disjoint_multi(bitset128!{'x', 'y'},       bitset128!{'z'},            bitset128!{'x', 'y', 'z'})]
  #[trace]
  fn test_bitset128_or_ref_right(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a | &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::both_empty(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'x', 'q'})]
  #[case::overlap_multi(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case::disjoint_multi(bitset128!{'x', 'y'},       bitset128!{'z'},            bitset128!{'x', 'y', 'z'})]
  #[trace]
  fn test_bitset128_or_ref_left(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a | b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::both_empty(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'x', 'q'})]
  #[case::overlap_multi(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case::disjoint_multi(bitset128!{'x', 'y'},       bitset128!{'z'},            bitset128!{'x', 'y', 'z'})]
  #[trace]
  fn test_bitset128_or_ref_both(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a | &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::both_empty(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'x', 'q'})]
  #[case::overlap_multi(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case::disjoint_multi(bitset128!{'x', 'y'},       bitset128!{'z'},            bitset128!{'x', 'y', 'z'})]
  #[trace]
  fn test_bitset128_or_assign(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a |= b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case::both_empty(bitset128!{},               bitset128!{},               bitset128!{})]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a', 'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'x', 'q'})]
  #[case::overlap_multi(bitset128!{'m', 'n'},       bitset128!{'a', 'b', 'm'},  bitset128!{'a', 'b', 'm', 'n'})]
  #[case::disjoint_multi(bitset128!{'x', 'y'},       bitset128!{'z'},            bitset128!{'x', 'y', 'z'})]
  #[trace]
  fn test_bitset128_or_assign_ref(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a |= &b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case::add_new(bitset128!{'a'},            b'b',                        bitset128!{'a', 'b'})]
  #[case::add_multi(bitset128!{'x', 'y'},       b'z',                        bitset128!{'x', 'y', 'z'})]
  #[case::add_existing(bitset128!{'m', 'n'},       b'm',                        bitset128!{'m', 'n'})]
  #[case::add_to_empty(bitset128!{},               b'a',                        bitset128!{'a'})]
  #[trace]
  fn test_bitset128_or_char_to_set(#[case] a: BitSet128, #[case] b: u8, #[case] expected: BitSet128) {
    let actual = a | b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::add_new(b'a',                        bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::add_multi(b'x',                        bitset128!{'y', 'z'},       bitset128!{'x', 'y', 'z'})]
  #[case::add_existing(b'm',                        bitset128!{'m', 'n'},       bitset128!{'m', 'n'})]
  #[case::add_to_empty(b'a',                        bitset128!{},               bitset128!{'a'})]
  #[trace]
  fn test_bitset128_or_set_to_char(#[case] a: u8, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = b | a;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::left_empty(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{})]
  #[case::partial(bitset128!{'a', 'b'},       bitset128!{'b'},            bitset128!{'a'})]
  #[case::disjoint(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e'},       bitset128!{'a', 'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a'})]
  #[case::overlap_one(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'x'})]
  #[trace]
  fn test_bitset128_sub(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a - b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::left_empty(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{})]
  #[case::partial(bitset128!{'a', 'b'},       bitset128!{'b'},            bitset128!{'a'})]
  #[case::disjoint(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e'},       bitset128!{'a', 'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a'})]
  #[case::overlap_one(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'x'})]
  #[trace]
  fn test_bitset128_sub_ref_right(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a - &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::left_empty(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{})]
  #[case::partial(bitset128!{'a', 'b'},       bitset128!{'b'},            bitset128!{'a'})]
  #[case::disjoint(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e'},       bitset128!{'a', 'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a'})]
  #[case::overlap_one(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'x'})]
  #[trace]
  fn test_bitset128_sub_ref_left(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a - b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::left_empty(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{})]
  #[case::partial(bitset128!{'a', 'b'},       bitset128!{'b'},            bitset128!{'a'})]
  #[case::disjoint(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e'},       bitset128!{'a', 'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a'})]
  #[case::overlap_one(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'x'})]
  #[trace]
  fn test_bitset128_sub_ref_both(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a - &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::left_empty(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{})]
  #[case::partial(bitset128!{'a', 'b'},       bitset128!{'b'},            bitset128!{'a'})]
  #[case::disjoint(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e'},       bitset128!{'a', 'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a'})]
  #[case::overlap_one(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'x'})]
  #[trace]
  fn test_bitset128_sub_assign(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a -= b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::left_empty(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{})]
  #[case::partial(bitset128!{'a', 'b'},       bitset128!{'b'},            bitset128!{'a'})]
  #[case::disjoint(bitset128!{'a', 'b', 'c'},  bitset128!{'d', 'e'},       bitset128!{'a', 'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a'})]
  #[case::overlap_one(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'x'})]
  #[trace]
  fn test_bitset128_sub_assign_ref(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a -= &b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case::empty_set(bitset128!{},               b'b',                        bitset128!{})]
  #[case::not_present(bitset128!{'a'},            b'b',                        bitset128!{'a'})]
  #[case::remove_one(bitset128!{'x', 'y'},       b'y',                        bitset128!{'x'})]
  #[case::remove_first(bitset128!{'m', 'n'},       b'm',                        bitset128!{'n'})]
  #[case::remove_only(bitset128!{'a'},            b'a',                        bitset128!{})]
  #[trace]
  fn test_bitset128_sub_char_from_set(#[case] a: BitSet128, #[case] b: u8, #[case] expected: BitSet128) {
    let actual = a - b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_empty_right(bitset128!{'a'},            bitset128!{},               bitset128!{})]
  #[case::right_empty_left(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'x'})]
  #[case::partial_match(bitset128!{'a', 'b', 'c'},  bitset128!{'a', 'd', 'e'},  bitset128!{'a'})]
  #[case::overlap_multi(bitset128!{'m', 'n', 'o'},  bitset128!{'n', 'o', 'p'},  bitset128!{'n', 'o'})]
  #[case::overlap_single(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'y'})]
  #[case::no_overlap(bitset128!{'g', 'h', 'i'},  bitset128!{'a', 'b', 'c'},  bitset128!{})]
  #[trace]
  fn test_bitset128_and(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a & b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_empty_right(bitset128!{'a'},            bitset128!{},               bitset128!{})]
  #[case::right_empty_left(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'x'})]
  #[case::partial_match(bitset128!{'a', 'b', 'c'},  bitset128!{'a', 'd', 'e'},  bitset128!{'a'})]
  #[case::overlap_multi(bitset128!{'m', 'n', 'o'},  bitset128!{'n', 'o', 'p'},  bitset128!{'n', 'o'})]
  #[case::overlap_single(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'y'})]
  #[case::no_overlap(bitset128!{'g', 'h', 'i'},  bitset128!{'a', 'b', 'c'},  bitset128!{})]
  #[trace]
  fn test_bitset128_and_ref_right(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a & &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_empty_right(bitset128!{'a'},            bitset128!{},               bitset128!{})]
  #[case::right_empty_left(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'x'})]
  #[case::partial_match(bitset128!{'a', 'b', 'c'},  bitset128!{'a', 'd', 'e'},  bitset128!{'a'})]
  #[case::overlap_multi(bitset128!{'m', 'n', 'o'},  bitset128!{'n', 'o', 'p'},  bitset128!{'n', 'o'})]
  #[case::overlap_single(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'y'})]
  #[case::no_overlap(bitset128!{'g', 'h', 'i'},  bitset128!{'a', 'b', 'c'},  bitset128!{})]
  #[trace]
  fn test_bitset128_and_ref_left(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a & b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_empty_right(bitset128!{'a'},            bitset128!{},               bitset128!{})]
  #[case::right_empty_left(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'x'})]
  #[case::partial_match(bitset128!{'a', 'b', 'c'},  bitset128!{'a', 'd', 'e'},  bitset128!{'a'})]
  #[case::overlap_multi(bitset128!{'m', 'n', 'o'},  bitset128!{'n', 'o', 'p'},  bitset128!{'n', 'o'})]
  #[case::overlap_single(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'y'})]
  #[case::no_overlap(bitset128!{'g', 'h', 'i'},  bitset128!{'a', 'b', 'c'},  bitset128!{})]
  #[trace]
  fn test_bitset128_and_ref_both(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a & &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_empty_right(bitset128!{'a'},            bitset128!{},               bitset128!{})]
  #[case::right_empty_left(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'x'})]
  #[case::partial_match(bitset128!{'a', 'b', 'c'},  bitset128!{'a', 'd', 'e'},  bitset128!{'a'})]
  #[case::overlap_multi(bitset128!{'m', 'n', 'o'},  bitset128!{'n', 'o', 'p'},  bitset128!{'n', 'o'})]
  #[case::overlap_single(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'y'})]
  #[case::no_overlap(bitset128!{'g', 'h', 'i'},  bitset128!{'a', 'b', 'c'},  bitset128!{})]
  #[trace]
  fn test_bitset128_and_assign(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a &= b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case::left_empty_right(bitset128!{'a'},            bitset128!{},               bitset128!{})]
  #[case::right_empty_left(bitset128!{},               bitset128!{'b'},            bitset128!{})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{})]
  #[case::identical(bitset128!{'a'},            bitset128!{'a'},            bitset128!{'a'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'b', 'c'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'x'})]
  #[case::partial_match(bitset128!{'a', 'b', 'c'},  bitset128!{'a', 'd', 'e'},  bitset128!{'a'})]
  #[case::overlap_multi(bitset128!{'m', 'n', 'o'},  bitset128!{'n', 'o', 'p'},  bitset128!{'n', 'o'})]
  #[case::overlap_single(bitset128!{'x', 'y'},       bitset128!{'y', 'z'},       bitset128!{'y'})]
  #[case::no_overlap(bitset128!{'g', 'h', 'i'},  bitset128!{'a', 'b', 'c'},  bitset128!{})]
  #[trace]
  fn test_bitset128_and_assign_ref(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a &= &b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case::present(bitset128!{'a', 'b'},       b'a',                        bitset128!{'a'})]
  #[case::not_present(bitset128!{'x', 'y'},       b'z',                        bitset128!{})]
  #[case::present_multi(bitset128!{'m', 'n', 'o'},  b'o',                        bitset128!{'o'})]
  #[trace]
  fn test_bitset128_and_char_to_set(#[case] a: BitSet128, #[case] b: u8, #[case] expected: BitSet128) {
    let actual = a & b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::present(b'a',                        bitset128!{'b', 'a'},       bitset128!{'a'})]
  #[case::not_present(b'z',                        bitset128!{'x', 'y'},       bitset128!{})]
  #[case::present_multi(b'o',                        bitset128!{'m', 'n', 'o'},  bitset128!{'o'})]
  #[trace]
  fn test_bitset128_and_set_to_char(#[case] a: u8, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = b & a;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q'})]
  #[case::identical(bitset128!{'a', 'b'},       bitset128!{'a', 'b'},       bitset128!{})]
  #[case::overlap_single(bitset128!{'m', 'n'},       bitset128!{'n', 'o'},       bitset128!{'m', 'o'})]
  #[case::overlap_multi(bitset128!{'x', 'y', 'z'},  bitset128!{'y', 'z', 'a'},  bitset128!{'x', 'a'})]
  #[case::alternating(bitset128!{'a', 'c', 'e'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a', 'b', 'd', 'e'})]
  #[trace]
  fn test_bitset128_xor(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a ^ b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q'})]
  #[case::identical(bitset128!{'a', 'b'},       bitset128!{'a', 'b'},       bitset128!{})]
  #[case::overlap_single(bitset128!{'m', 'n'},       bitset128!{'n', 'o'},       bitset128!{'m', 'o'})]
  #[case::overlap_multi(bitset128!{'x', 'y', 'z'},  bitset128!{'y', 'z', 'a'},  bitset128!{'x', 'a'})]
  #[case::alternating(bitset128!{'a', 'c', 'e'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a', 'b', 'd', 'e'})]
  #[trace]
  fn test_bitset128_xor_ref_right(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = a ^ &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q'})]
  #[case::identical(bitset128!{'a', 'b'},       bitset128!{'a', 'b'},       bitset128!{})]
  #[case::overlap_single(bitset128!{'m', 'n'},       bitset128!{'n', 'o'},       bitset128!{'m', 'o'})]
  #[case::overlap_multi(bitset128!{'x', 'y', 'z'},  bitset128!{'y', 'z', 'a'},  bitset128!{'x', 'a'})]
  #[case::alternating(bitset128!{'a', 'c', 'e'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a', 'b', 'd', 'e'})]
  #[trace]
  fn test_bitset128_xor_ref_left(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a ^ b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q'})]
  #[case::identical(bitset128!{'a', 'b'},       bitset128!{'a', 'b'},       bitset128!{})]
  #[case::overlap_single(bitset128!{'m', 'n'},       bitset128!{'n', 'o'},       bitset128!{'m', 'o'})]
  #[case::overlap_multi(bitset128!{'x', 'y', 'z'},  bitset128!{'y', 'z', 'a'},  bitset128!{'x', 'a'})]
  #[case::alternating(bitset128!{'a', 'c', 'e'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a', 'b', 'd', 'e'})]
  #[trace]
  fn test_bitset128_xor_ref_both(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = &a ^ &b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q'})]
  #[case::identical(bitset128!{'a', 'b'},       bitset128!{'a', 'b'},       bitset128!{})]
  #[case::overlap_single(bitset128!{'m', 'n'},       bitset128!{'n', 'o'},       bitset128!{'m', 'o'})]
  #[case::overlap_multi(bitset128!{'x', 'y', 'z'},  bitset128!{'y', 'z', 'a'},  bitset128!{'x', 'a'})]
  #[case::alternating(bitset128!{'a', 'c', 'e'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a', 'b', 'd', 'e'})]
  #[trace]
  fn test_bitset128_xor_assign(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a ^= b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case::left_only(bitset128!{'a'},            bitset128!{},               bitset128!{'a'})]
  #[case::right_only(bitset128!{},               bitset128!{'b'},            bitset128!{'b'})]
  #[case::disjoint(bitset128!{'a'},            bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::superset(bitset128!{'a', 'b', 'c'},  bitset128!{'b', 'c'},       bitset128!{'a'})]
  #[case::partial_overlap(bitset128!{'a', 'b', 'x'},  bitset128!{'p', 'x', 'q'},  bitset128!{'a', 'b', 'p', 'q'})]
  #[case::identical(bitset128!{'a', 'b'},       bitset128!{'a', 'b'},       bitset128!{})]
  #[case::overlap_single(bitset128!{'m', 'n'},       bitset128!{'n', 'o'},       bitset128!{'m', 'o'})]
  #[case::overlap_multi(bitset128!{'x', 'y', 'z'},  bitset128!{'y', 'z', 'a'},  bitset128!{'x', 'a'})]
  #[case::alternating(bitset128!{'a', 'c', 'e'},  bitset128!{'b', 'c', 'd'},  bitset128!{'a', 'b', 'd', 'e'})]
  #[trace]
  fn test_bitset128_xor_assign_ref(#[case] a: BitSet128, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let mut a = a;
    a ^= &b;
    assert_eq!(a, expected);
  }

  #[rstest]
  #[case::add_new(bitset128!{'a'},            b'b',                        bitset128!{'a', 'b'})]
  #[case::toggle_present(bitset128!{'x', 'y'},       b'x',                        bitset128!{'y'})]
  #[case::toggle_first(bitset128!{'m', 'n'},       b'm',                        bitset128!{'n'})]
  #[case::toggle_second(bitset128!{'c', 'd'},       b'd',                        bitset128!{'c'})]
  #[trace]
  fn test_bitset128_xor_char_to_set(#[case] a: BitSet128, #[case] b: u8, #[case] expected: BitSet128) {
    let actual = a ^ b;
    assert_eq!(actual, expected);
  }

  #[rstest]
  #[case::add_new(b'a',                        bitset128!{'b'},            bitset128!{'a', 'b'})]
  #[case::toggle_present(b'x',                        bitset128!{'x', 'y'},       bitset128!{'y'})]
  #[case::toggle_first(b'm',                        bitset128!{'m', 'n'},       bitset128!{'n'})]
  #[case::toggle_second(b'd',                        bitset128!{'c', 'd'},       bitset128!{'c'})]
  #[trace]
  fn test_bitset128_xor_set_to_char(#[case] a: u8, #[case] b: BitSet128, #[case] expected: BitSet128) {
    let actual = b ^ a;
    assert_eq!(actual, expected);
  }
}
