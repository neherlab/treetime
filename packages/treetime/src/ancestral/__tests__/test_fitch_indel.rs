#[cfg(test)]
mod tests {
  use crate::ancestral::fitch_indel::{resolve_indels_backward, resolve_indels_forward};
  use pretty_assertions::assert_eq;
  use std::collections::BTreeSet;
  use treetime_primitives::Seq;

  fn refs<T>(v: &[T]) -> Vec<&T> {
    v.iter().collect()
  }

  #[test]
  fn test_fitch_indel_backward_no_disagreement() {
    let child_gaps = [vec![(2, 4)], vec![(2, 4)]];
    let child_unknown = [vec![], vec![]];
    let empty0 = BTreeSet::new();
    let empty1 = BTreeSet::new();
    let child_vis = [&empty0, &empty1];
    let result = resolve_indels_backward(&refs(&child_gaps), &refs(&child_unknown), &child_vis, 10);
    assert!(
      result.variable_indel.is_empty(),
      "No disagreement when all children have same gaps"
    );
    assert_eq!(result.resolved_gaps, vec![(2, 4)]);
  }

  #[test]
  fn test_fitch_indel_backward_one_child_has_gap() {
    let child_gaps = [vec![(2, 4)], vec![]];
    let child_unknown = [vec![], vec![]];
    let empty0 = BTreeSet::new();
    let empty1 = BTreeSet::new();
    let child_vis = [&empty0, &empty1];
    let result = resolve_indels_backward(&refs(&child_gaps), &refs(&child_unknown), &child_vis, 10);
    assert_eq!(result.variable_indel.len(), 1);
    assert!(result.variable_indel.contains(&(2, 4)));
  }

  #[test]
  fn test_fitch_indel_backward_all_children_deleted() {
    let child_gaps = [vec![(2, 4)], vec![(2, 4)], vec![(2, 4)]];
    let child_unknown = [vec![], vec![], vec![]];
    let empty0 = BTreeSet::new();
    let empty1 = BTreeSet::new();
    let empty2 = BTreeSet::new();
    let child_vis = [&empty0, &empty1, &empty2];
    let result = resolve_indels_backward(&refs(&child_gaps), &refs(&child_unknown), &child_vis, 10);
    assert!(
      result.variable_indel.is_empty(),
      "All children agree on gap, no variable indels"
    );
    assert_eq!(result.resolved_gaps, vec![(2, 4)]);
  }

  #[test]
  fn test_fitch_indel_backward_propagates_child_variable_indels() {
    let child_gaps = [vec![(2, 4)], vec![]];
    let child_unknown = [vec![], vec![]];
    let empty0 = BTreeSet::new();
    let child1_vi = BTreeSet::from([(2_usize, 4_usize)]);
    let child_vis = [&empty0, &child1_vi];
    let result = resolve_indels_backward(&refs(&child_gaps), &refs(&child_unknown), &child_vis, 10);

    // Child 0 has gap, child 1 has variable_indel: both compatible with gap -> resolved.
    assert!(result.variable_indel.is_empty());
    assert_eq!(result.resolved_gaps, vec![(2, 4)]);
  }

  #[test]
  fn test_fitch_indel_backward_all_unknown_not_deletion_evidence() {
    let child_gaps = [vec![], vec![]];
    let child_unknown = [vec![(2, 4)], vec![(2, 4)]];
    let empty0 = BTreeSet::new();
    let empty1 = BTreeSet::new();
    let child_vis = [&empty0, &empty1];
    let result = resolve_indels_backward(&refs(&child_gaps), &refs(&child_unknown), &child_vis, 10);
    assert!(result.variable_indel.is_empty());
    assert!(result.resolved_gaps.is_empty());
  }

  #[test]
  fn test_fitch_indel_backward_propagates_child_variable_indel_set() {
    let child_gaps = [vec![], vec![]];
    let child_unknown = [vec![], vec![]];
    let vi0 = BTreeSet::from([(2_usize, 4_usize)]);
    let vi1 = BTreeSet::from([(2_usize, 4_usize)]);
    let child_vis = [&vi0, &vi1];
    let result = resolve_indels_backward(&refs(&child_gaps), &refs(&child_unknown), &child_vis, 10);
    assert_eq!(result.variable_indel.len(), 1);
    assert!(result.resolved_gaps.is_empty());
  }

  #[test]
  fn test_fitch_indel_backward_single_child_preserves_gaps() {
    let child_gaps = [vec![(2, 4)]];
    let child_unknown = [vec![]];
    let empty0 = BTreeSet::new();
    let child_vis = [&empty0];
    let result = resolve_indels_backward(&refs(&child_gaps), &refs(&child_unknown), &child_vis, 10);
    assert!(result.variable_indel.is_empty());
    assert_eq!(result.resolved_gaps, vec![(2, 4)]);
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
    assert!(indels[0].deletion);
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
    assert!(!indels[0].deletion);
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
