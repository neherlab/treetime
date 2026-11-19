#[cfg(test)]
mod tests {
  use crate::seq::indel::{resolve_indels_backward, resolve_indels_forward};
  use itertools::Itertools;
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use treetime_primitives::Seq;

  #[rustfmt::skip]
  #[rstest]
  #[case::children_agree_on_gap(         vec![vec![(2, 4)], vec![(2, 4)]], vec![vec![],       vec![]],       vec![btreeset! {},       btreeset! {}],       (vec![(2, 4)], btreeset! {}))]
  #[case::one_child_has_gap(             vec![vec![(2, 4)], vec![]],       vec![vec![],       vec![]],       vec![btreeset! {},       btreeset! {}],       (vec![],       btreeset! {(2, 4)}))]
  #[case::gap_meets_child_variable_indel(vec![vec![(2, 4)], vec![]],       vec![vec![],       vec![]],       vec![btreeset! {},       btreeset! {(2, 4)}], (vec![(2, 4)], btreeset! {}))]
  #[case::all_children_unknown(          vec![vec![],       vec![]],       vec![vec![(2, 4)], vec![(2, 4)]], vec![btreeset! {},       btreeset! {}],       (vec![],       btreeset! {}))]
  #[case::all_children_variable(         vec![vec![],       vec![]],       vec![vec![],       vec![]],       vec![btreeset! {(2, 4)}, btreeset! {(2, 4)}], (vec![],       btreeset! {(2, 4)}))]
  #[case::single_child_keeps_its_gaps(   vec![vec![(2, 4)]],               vec![vec![]],                     vec![btreeset! {}],                         (vec![(2, 4)], btreeset! {}))]
  #[trace]
  fn test_fitch_indel_backward_resolves_gaps_and_variable_indels(
    #[case] child_gaps: Vec<Vec<(usize, usize)>>,
    #[case] child_unknown: Vec<Vec<(usize, usize)>>,
    #[case] child_variable_indels: Vec<BTreeSet<(usize, usize)>>,
    #[case] expected: (Vec<(usize, usize)>, BTreeSet<(usize, usize)>),
  ) {
    let result = resolve_indels_backward(
      &child_gaps.iter().collect_vec(),
      &child_unknown.iter().collect_vec(),
      &child_variable_indels.iter().collect_vec(),
      10,
    );
    assert_eq!(expected, (result.resolved_gaps, result.variable_indel));
  }

  #[test]
  fn test_fitch_indel_forward_variable_parent_has_sequence() {
    let variable_indel = BTreeSet::from([(2_usize, 5_usize)]);

    let parent_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let node_seq = Seq::try_from_str("ACGTACGTAC").unwrap();

    let (indels, new_gaps) = resolve_indels_forward(&variable_indel, &[], &[], &[], &parent_seq, &node_seq);

    assert!(indels.is_empty());
    assert!(new_gaps.is_empty());
  }

  #[test]
  fn test_fitch_indel_forward_variable_parent_has_gap() {
    let variable_indel = BTreeSet::from([(2_usize, 5_usize)]);

    let parent_gaps = [(2, 5)];
    let parent_seq = Seq::try_from_str("AC---CGTAC").unwrap();
    let node_seq = Seq::try_from_str("ACGTACGTAC").unwrap();

    let (indels, new_gaps) = resolve_indels_forward(&variable_indel, &[], &[], &parent_gaps, &parent_seq, &node_seq);

    assert!(indels.is_empty());
    assert_eq!(new_gaps, vec![(2, 5)]);
  }

  #[test]
  fn test_fitch_indel_forward_consensus_deletion() {
    let variable_indel = BTreeSet::new();
    let parent_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let node_seq = Seq::try_from_str("AC---CGTAC").unwrap();
    let node_gaps = [(2, 5)];

    let (indels, new_gaps) = resolve_indels_forward(&variable_indel, &node_gaps, &[], &[], &parent_seq, &node_seq);

    assert_eq!(indels.len(), 1);
    assert!(indels[0].is_deletion());
    assert_eq!(indels[0].range, (2, 5));
    assert_eq!(indels[0].seq, Seq::try_from_str("GTA").unwrap());
    assert_eq!(new_gaps, vec![(2, 5)]);
  }

  #[test]
  fn test_fitch_indel_forward_consensus_insertion() {
    let variable_indel = BTreeSet::new();
    let parent_gaps = [(2, 5)];
    let parent_seq = Seq::try_from_str("AC---CGTAC").unwrap();
    let node_seq = Seq::try_from_str("ACGTACGTAC").unwrap();

    let (indels, new_gaps) = resolve_indels_forward(&variable_indel, &[], &[], &parent_gaps, &parent_seq, &node_seq);

    assert_eq!(indels.len(), 1);
    assert!(!indels[0].is_deletion());
    assert_eq!(indels[0].range, (2, 5));
    assert_eq!(indels[0].seq, Seq::try_from_str("GTA").unwrap());
    assert!(new_gaps.is_empty());
  }

  #[test]
  fn test_fitch_indel_forward_no_indels() {
    let variable_indel = BTreeSet::new();
    let parent_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let node_seq = Seq::try_from_str("ACGTACGTAC").unwrap();

    let (indels, _) = resolve_indels_forward(&variable_indel, &[], &[], &[], &parent_seq, &node_seq);

    assert!(indels.is_empty());
  }
}
