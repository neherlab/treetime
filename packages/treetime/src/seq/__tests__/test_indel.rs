#[cfg(test)]
mod tests {
  use crate::seq::indel::{InDel, compose_indels, sort_indels};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_primitives::Seq;

  fn del(start: usize, end: usize, seq: &str) -> InDel {
    InDel::del((start, end), Seq::try_from_str(seq).unwrap()).unwrap()
  }

  fn ins(start: usize, end: usize, seq: &str) -> InDel {
    InDel::ins((start, end), Seq::try_from_str(seq).unwrap()).unwrap()
  }

  #[test]
  fn test_indel_compose_both_empty() {
    let result = compose_indels(&[], &[]);
    assert_eq!(result, vec![]);
  }

  #[test]
  fn test_indel_compose_parent_empty() {
    let child = vec![del(1, 4, "CGT")];
    let result = compose_indels(&[], &child);
    assert_eq!(result, child);
  }

  #[test]
  fn test_indel_compose_child_empty() {
    let parent = vec![del(1, 4, "CGT")];
    let result = compose_indels(&parent, &[]);
    assert_eq!(result, parent);
  }

  #[test]
  fn test_indel_compose_case1_overlapping_deletions_same_range() {
    let parent = vec![del(1, 3, "CG")];
    let child = vec![del(1, 4, "CGT")];
    let result = compose_indels(&parent, &child);
    let expected = vec![del(1, 4, "CGT")];
    assert_eq!(result, expected);
  }

  #[test]
  fn test_indel_compose_case2_adjacent_deletions() {
    let parent = vec![del(1, 3, "CG")];
    let child = vec![del(3, 5, "TC")];
    let result = compose_indels(&parent, &child);
    let expected = vec![del(1, 5, "CGTC")];
    assert_eq!(result, expected);
  }

  #[test]
  fn test_indel_compose_case3_overlapping_insertions_same_range() {
    let parent = vec![ins(1, 4, "TTT")];
    let child = vec![ins(1, 4, "GGG")];
    let result = compose_indels(&parent, &child);
    let expected = vec![ins(1, 4, "GGG")];
    assert_eq!(result, expected);
  }

  #[test]
  fn test_indel_compose_case4_deletion_then_insertion_cancel() {
    let parent = vec![del(1, 4, "CGT")];
    let child = vec![ins(1, 4, "CGT")];
    let result = compose_indels(&parent, &child);
    let expected: Vec<InDel> = vec![];
    assert_eq!(result, expected);
  }

  #[test]
  fn test_indel_compose_case5_deletion_then_insertion_different() {
    let parent = vec![del(1, 4, "CGT")];
    let child = vec![ins(1, 4, "TTT")];
    let result = compose_indels(&parent, &child);
    let expected = vec![del(1, 4, "CGT"), ins(1, 4, "TTT")];
    assert_eq!(result, expected);
  }

  #[test]
  fn test_indel_compose_case6_partially_overlapping_deletions() {
    // Overlap at 2..3. Parent provides seq for overlap, child for non-overlap suffix.
    // Seq assembled from available data, not ground-truth ancestor content.
    let parent = vec![del(1, 3, "CG")];
    let child = vec![del(2, 5, "TAC")];
    let result = compose_indels(&parent, &child);
    let expected = vec![del(1, 5, "CGAC")];
    assert_eq!(result, expected);
  }

  #[test]
  fn test_indel_compose_case7_non_overlapping() {
    let parent = vec![del(1, 3, "CG")];
    let child = vec![del(6, 7, "G")];
    let result = compose_indels(&parent, &child);
    let expected = vec![del(1, 3, "CG"), del(6, 7, "G")];
    assert_eq!(result, expected);
  }

  #[test]
  fn test_indel_compose_multiple_non_overlapping_interleaved() {
    let parent = vec![del(1, 3, "CG"), del(10, 12, "AT")];
    let child = vec![del(5, 7, "GC"), del(15, 16, "T")];
    let result = compose_indels(&parent, &child);
    let expected = vec![del(1, 3, "CG"), del(5, 7, "GC"), del(10, 12, "AT"), del(15, 16, "T")];
    assert_eq!(result, expected);
  }

  #[test]
  fn test_indel_compose_insertion_then_deletion_cancel() {
    let parent = vec![ins(2, 5, "AAA")];
    let child = vec![del(2, 5, "AAA")];
    let result = compose_indels(&parent, &child);
    let expected: Vec<InDel> = vec![];
    assert_eq!(result, expected);
  }

  #[test]
  fn test_indel_compose_insertion_then_deletion_different() {
    let parent = vec![ins(2, 5, "AAA")];
    let child = vec![del(2, 5, "GGG")];
    let result = compose_indels(&parent, &child);
    let expected = vec![ins(2, 5, "AAA"), del(2, 5, "GGG")];
    assert_eq!(result, expected);
  }

  #[test]
  fn test_indel_compose_mixed_overlap_output_sorted() {
    // del at higher range, ins at lower range: output must be sorted
    let parent = vec![del(3, 7, "ACGT")];
    let child = vec![ins(1, 5, "TTTT")];
    let result = compose_indels(&parent, &child);
    assert_eq!(result, vec![ins(1, 5, "TTTT"), del(3, 7, "ACGT")]);
    assert!(
      result.windows(2).all(|w| match w {
        [a, b] => a.range.0 <= b.range.0,
        _ => true,
      }),
      "output must be sorted by range.0"
    );
  }

  #[test]
  fn test_indel_compose_child_before_parent_adjacent() {
    // Child deletion ends where parent starts - merged by final pass
    let parent = vec![del(5, 8, "ACG")];
    let child = vec![del(3, 5, "TT")];
    let result = compose_indels(&parent, &child);
    let expected = vec![del(3, 8, "TTACG")];
    assert_eq!(result, expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::adjacent_triple(
    vec![del(0, 2, "AC")],
    vec![del(2, 4, "GT"), del(4, 6, "CG")],
    vec![del(0, 6, "ACGTCG")]
  )]
  #[case::parent_superset(
    vec![del(0, 10, "ACGTACGTAC")],
    vec![del(2, 5, "TAC")],
    vec![del(0, 10, "ACGTACGTAC")]
  )]
  #[trace]
  fn test_indel_compose_parametric(
    #[case] parent: Vec<InDel>,
    #[case] child: Vec<InDel>,
    #[case] expected: Vec<InDel>,
  ) {
    let result = compose_indels(&parent, &child);
    assert_eq!(result, expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::del_after_ins(  vec![del(5, 8, "ACG")],                vec![ins(1, 4, "TTT")])]
  #[case::ins_before_del( vec![ins(1, 4, "TTT")],                vec![del(5, 8, "ACG")])]
  #[case::overlap_del_ins(vec![del(1, 5, "ACGT")],               vec![ins(3, 7, "TTTT")])]
  #[case::overlap_ins_del(vec![ins(3, 7, "TTTT")],               vec![del(1, 5, "ACGT")])]
  #[case::multi_del(      vec![del(1, 3, "CG"), del(6, 8, "AT")], vec![del(4, 5, "T")])]
  #[trace]
  fn test_indel_compose_output_sorted_invariant(#[case] parent: Vec<InDel>, #[case] child: Vec<InDel>) {
    let result = compose_indels(&parent, &child);
    for w in result.windows(2) {
      let [a, b] = w else { continue };
      assert!(
        a.range.0 <= b.range.0,
        "output not sorted: result={result:?}"
      );
    }
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::short(  (0,   3),   "ACG")]
  #[case::medium( (5,   10),  "ACGTC")]
  #[case::long(   (100, 105), "TTTTT")]
  #[trace]
  fn test_indel_compose_cancellation_roundtrip(#[case] (start, end): (usize, usize), #[case] seq: &str) {
    let result = compose_indels(&[del(start, end, seq)], &[ins(start, end, seq)]);
    assert_eq!(result, vec![], "del+ins should cancel");

    let result = compose_indels(&[ins(start, end, seq)], &[del(start, end, seq)]);
    assert_eq!(result, vec![], "ins+del should cancel");
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::partial_overlap(  vec![del(1, 3, "CG")],  vec![del(2, 5, "TAC")])]
  #[case::triple_merge(     vec![del(0, 2, "AC")],   vec![del(2, 4, "GT"), del(4, 6, "CG")])]
  #[case::del_then_ins(     vec![del(1, 4, "CGT")],  vec![ins(1, 4, "TTT")])]
  #[case::ins_then_ins(     vec![ins(1, 4, "TTT")],  vec![ins(1, 4, "GGG")])]
  #[trace]
  fn test_indel_compose_seq_length_invariant(#[case] parent: Vec<InDel>, #[case] child: Vec<InDel>) {
    let result = compose_indels(&parent, &child);
    for indel in &result {
      assert_eq!(
        indel.seq.len(),
        indel.range.1 - indel.range.0,
        "seq length mismatch for {indel:?}"
      );
    }
  }

  #[test]
  fn test_indel_sort_indels() {
    let mut indels = vec![del(5, 7, "AC"), del(1, 3, "CG"), del(3, 5, "GT")];
    sort_indels(&mut indels);
    assert_eq!(indels, vec![del(1, 3, "CG"), del(3, 5, "GT"), del(5, 7, "AC")]);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::reversed((3, 2), "A", "non-empty and ordered")]
  #[case::empty(   (2, 2), "",  "non-empty and ordered")]
  #[case::short(   (2, 4), "A", "sequence has length 1")]
  #[trace]
  fn test_indel_constructor_rejects_invalid_event(
    #[case] range: (usize, usize),
    #[case] sequence: &str,
    #[case] message: &str,
  ) {
    let sequence = Seq::try_from_str(sequence).unwrap();
    let error = InDel::del(range, sequence).expect_err("invalid indel must be rejected");
    assert!(error.to_string().contains(message));
  }
}
