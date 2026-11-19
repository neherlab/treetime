#[cfg(test)]
mod tests {
  use crate::seq::indel::{InDel, compose_indels, sort_indels};
  use helpers::{del, ins};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::both_empty(                     vec![],                                  vec![],                                  vec![])]
  #[case::parent_empty(                   vec![],                                  vec![del(1, 4, "CGT")],                  vec![del(1, 4, "CGT")])]
  #[case::child_empty(                    vec![del(1, 4, "CGT")],                  vec![],                                  vec![del(1, 4, "CGT")])]
  #[case::deletions_on_the_same_start(    vec![del(1, 3, "CG")],                   vec![del(1, 4, "CGT")],                  vec![del(1, 4, "CGT")])]
  #[case::adjacent_deletions(             vec![del(1, 3, "CG")],                   vec![del(3, 5, "TC")],                   vec![del(1, 5, "CGTC")])]
  #[case::child_deletion_before_parent(   vec![del(5, 8, "ACG")],                  vec![del(3, 5, "TT")],                   vec![del(3, 8, "TTACG")])]
  #[case::adjacent_triple(                vec![del(0, 2, "AC")],                   vec![del(2, 4, "GT"), del(4, 6, "CG")],  vec![del(0, 6, "ACGTCG")])]
  #[case::partially_overlapping_deletions(vec![del(1, 3, "CG")],                   vec![del(2, 5, "TAC")],                  vec![del(1, 5, "CGAC")])]
  #[case::parent_deletion_contains_child( vec![del(0, 10, "ACGTACGTAC")],          vec![del(2, 5, "TAC")],                  vec![del(0, 10, "ACGTACGTAC")])]
  #[case::insertions_on_the_same_range(   vec![ins(1, 4, "TTT")],                  vec![ins(1, 4, "GGG")],                  vec![ins(1, 4, "GGG")])]
  #[case::deletion_then_same_insertion(   vec![del(1, 4, "CGT")],                  vec![ins(1, 4, "CGT")],                  vec![])]
  #[case::insertion_then_same_deletion(   vec![ins(2, 5, "AAA")],                  vec![del(2, 5, "AAA")],                  vec![])]
  #[case::deletion_then_other_insertion(  vec![del(1, 4, "CGT")],                  vec![ins(1, 4, "TTT")],                  vec![del(1, 4, "CGT"), ins(1, 4, "TTT")])]
  #[case::insertion_then_other_deletion(  vec![ins(2, 5, "AAA")],                  vec![del(2, 5, "GGG")],                  vec![ins(2, 5, "AAA"), del(2, 5, "GGG")])]
  #[case::overlapping_deletion_insertion( vec![del(3, 7, "ACGT")],                 vec![ins(1, 5, "TTTT")],                 vec![ins(1, 5, "TTTT"), del(3, 7, "ACGT")])]
  #[case::disjoint(                       vec![del(1, 3, "CG")],                   vec![del(6, 7, "G")],                    vec![del(1, 3, "CG"), del(6, 7, "G")])]
  #[case::disjoint_interleaved(           vec![del(1, 3, "CG"), del(10, 12, "AT")], vec![del(5, 7, "GC"), del(15, 16, "T")], vec![del(1, 3, "CG"), del(5, 7, "GC"), del(10, 12, "AT"), del(15, 16, "T")])]
  #[trace]
  fn test_indel_compose_applies_the_child_indels_after_the_parent_indels(
    #[case] parent: Vec<InDel>,
    #[case] child: Vec<InDel>,
    #[case] expected: Vec<InDel>,
  ) {
    let actual = compose_indels(&parent, &child);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::del_after_ins(  vec![del(5, 8, "ACG")],                vec![ins(1, 4, "TTT")])]
  #[case::ins_before_del( vec![ins(1, 4, "TTT")],                vec![del(5, 8, "ACG")])]
  #[case::overlap_del_ins(vec![del(1, 5, "ACGT")],               vec![ins(3, 7, "TTTT")])]
  #[case::overlap_ins_del(vec![ins(3, 7, "TTTT")],               vec![del(1, 5, "ACGT")])]
  #[case::multi_del(      vec![del(1, 3, "CG"), del(6, 8, "AT")], vec![del(4, 5, "T")])]
  #[trace]
  fn test_indel_compose_output_is_sorted_by_start_position(#[case] parent: Vec<InDel>, #[case] child: Vec<InDel>) {
    let result = compose_indels(&parent, &child);
    assert!(
      result.is_sorted_by_key(|x| x.range.0),
      "output not sorted: result={result:?}"
    );
  }

  #[test]
  fn test_indel_sort_indels_orders_by_start_position() {
    let mut actual = vec![del(5, 7, "AC"), del(1, 3, "CG"), del(3, 5, "GT")];
    sort_indels(&mut actual);
    assert_eq!(vec![del(1, 3, "CG"), del(3, 5, "GT"), del(5, 7, "AC")], actual);
  }

  mod helpers {
    use crate::seq::indel::InDel;
    use crate::test_utils::{deletion, insertion};
    use treetime_primitives::Seq;

    pub(super) fn del(start: usize, end: usize, seq: &str) -> InDel {
      deletion((start, end), Seq::try_from_str(seq).unwrap())
    }

    pub(super) fn ins(start: usize, end: usize, seq: &str) -> InDel {
      insertion((start, end), Seq::try_from_str(seq).unwrap())
    }
  }
}
