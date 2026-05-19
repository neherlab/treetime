#[cfg(test)]
mod tests {
  use crate::ancestral::fitch_indel::{resolve_indels_backward, resolve_indels_forward};
  use crate::representation::payload::sparse::Deletion;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_primitives::Seq;

  // -- H1 reproduction: all-unknown children must not become deletion evidence --

  #[test]
  fn test_fitch_indel_backward_all_unknown_not_deletion_evidence() {
    // Two children, both unknown at (2,4), no gaps anywhere.
    // Unknown = missing data, not observed deletion.
    // Expected: no resolved gaps, no variable indels.
    let child_gaps = vec![vec![], vec![]];
    let child_unknown = vec![vec![(2, 4)], vec![(2, 4)]];
    let empty0 = BTreeMap::new();
    let empty1 = BTreeMap::new();
    let child_vis: Vec<&BTreeMap<(usize, usize), Deletion>> = vec![&empty0, &empty1];
    let child_gaps_ref: Vec<&Vec<_>> = child_gaps.iter().collect();
    let child_unknown_ref: Vec<&Vec<_>> = child_unknown.iter().collect();

    let result = resolve_indels_backward(&child_gaps_ref, &child_unknown_ref, &child_vis, 10);

    assert!(result.resolved_gaps.is_empty());
    assert!(result.variable_indel.is_empty(), "All-unknown interval must not create deletion evidence");
  }

  // -- H2 reproduction: child variable-indel counts must propagate upward --

  #[test]
  fn test_fitch_indel_backward_propagates_variable_indel_counts() {
    // Two children, both have variable_indel at (2,4) with deleted=3, present=1.
    // No direct gaps or unknowns.
    // Expected: parent variable_indel at (2,4) aggregates child counts (deleted=6, present=2).
    let child_gaps = vec![vec![], vec![]];
    let child_unknown = vec![vec![], vec![]];

    let mut child0_vi = BTreeMap::new();
    child0_vi.insert((2, 4), Deletion { deleted: 3, present: 1 });

    let mut child1_vi = BTreeMap::new();
    child1_vi.insert((2, 4), Deletion { deleted: 3, present: 1 });

    let child_vis: Vec<&BTreeMap<(usize, usize), Deletion>> = vec![&child0_vi, &child1_vi];
    let child_gaps_ref: Vec<&Vec<_>> = child_gaps.iter().collect();
    let child_unknown_ref: Vec<&Vec<_>> = child_unknown.iter().collect();

    let result = resolve_indels_backward(&child_gaps_ref, &child_unknown_ref, &child_vis, 10);

    assert!(result.resolved_gaps.is_empty());
    let del = &result.variable_indel[&(2, 4)];
    assert!(del.deleted > del.present, "Deletion majority from child subtrees must be preserved, got deleted={}, present={}", del.deleted, del.present);
  }

  // -- H3 reproduction: single-child nodes must preserve fixed gap ranges --

  #[test]
  fn test_fitch_indel_backward_single_child_preserves_gaps() {
    // One child with a fixed gap at (2,4).
    // Expected: parent resolved_gaps = [(2,4)] (identity behavior for single child).
    let child_gaps = vec![vec![(2, 4)]];
    let child_unknown = vec![vec![]];
    let empty0 = BTreeMap::new();
    let child_vis: Vec<&BTreeMap<(usize, usize), Deletion>> = vec![&empty0];
    let child_gaps_ref: Vec<&Vec<_>> = child_gaps.iter().collect();
    let child_unknown_ref: Vec<&Vec<_>> = child_unknown.iter().collect();

    let result = resolve_indels_backward(&child_gaps_ref, &child_unknown_ref, &child_vis, 10);

    assert!(result.variable_indel.is_empty());
    assert_eq!(vec![(2, 4)], result.resolved_gaps, "Single-child gaps must propagate to parent");
  }

  #[test]
  fn test_fitch_indel_backward_no_disagreement() {
    let child_gaps = vec![vec![(2, 4)], vec![(2, 4)]];
    let child_unknown = vec![vec![], vec![]];
    let empty0 = BTreeMap::new();
    let empty1 = BTreeMap::new();
    let child_vis: Vec<&BTreeMap<(usize, usize), Deletion>> = vec![&empty0, &empty1];
    let child_gaps_ref: Vec<&Vec<_>> = child_gaps.iter().collect();
    let child_unknown_ref: Vec<&Vec<_>> = child_unknown.iter().collect();
    let result = resolve_indels_backward(&child_gaps_ref, &child_unknown_ref, &child_vis, 10);
    assert!(result.variable_indel.is_empty(), "No disagreement when all children have same gaps");
    assert_eq!(result.resolved_gaps, vec![(2, 4)]);
  }

  #[test]
  fn test_fitch_indel_backward_one_child_has_gap() {
    let child_gaps = vec![vec![(2, 4)], vec![]];
    let child_unknown = vec![vec![], vec![]];
    let empty0 = BTreeMap::new();
    let empty1 = BTreeMap::new();
    let child_vis: Vec<&BTreeMap<(usize, usize), Deletion>> = vec![&empty0, &empty1];
    let child_gaps_ref: Vec<&Vec<_>> = child_gaps.iter().collect();
    let child_unknown_ref: Vec<&Vec<_>> = child_unknown.iter().collect();
    let result = resolve_indels_backward(&child_gaps_ref, &child_unknown_ref, &child_vis, 10).variable_indel;
    assert_eq!(result.len(), 1);
    let del = &result[&(2, 4)];
    assert_eq!(del.deleted, 1);
    assert_eq!(del.present, 1);
  }

  #[test]
  fn test_fitch_indel_backward_all_children_deleted() {
    let child_gaps = vec![vec![(2, 4)], vec![(2, 4)], vec![(2, 4)]];
    let child_unknown = vec![vec![], vec![], vec![]];
    let empty0 = BTreeMap::new();
    let empty1 = BTreeMap::new();
    let empty2 = BTreeMap::new();
    let child_vis: Vec<&BTreeMap<(usize, usize), Deletion>> = vec![&empty0, &empty1, &empty2];
    let child_gaps_ref: Vec<&Vec<_>> = child_gaps.iter().collect();
    let child_unknown_ref: Vec<&Vec<_>> = child_unknown.iter().collect();

    let result = resolve_indels_backward(&child_gaps_ref, &child_unknown_ref, &child_vis, 10);
    assert!(result.variable_indel.is_empty(), "All children agree on gap, no variable indels");
    assert_eq!(result.resolved_gaps, vec![(2, 4)]);
  }

  #[test]
  fn test_fitch_indel_backward_propagates_child_variable_indels() {
    // Child 0 has gap at (2,4), child 1 does not.
    // Child 1 has variable indel at (2,4).
    let child_gaps = vec![vec![(2, 4)], vec![]];
    let child_unknown = vec![vec![], vec![]];
    let empty0 = BTreeMap::new();
    let mut child1_vi = BTreeMap::new();
    child1_vi.insert((2, 4), Deletion { deleted: 1, present: 1 });
    let child_vis: Vec<&BTreeMap<(usize, usize), Deletion>> = vec![&empty0, &child1_vi];
    let child_gaps_ref: Vec<&Vec<_>> = child_gaps.iter().collect();
    let child_unknown_ref: Vec<&Vec<_>> = child_unknown.iter().collect();
    let result = resolve_indels_backward(&child_gaps_ref, &child_unknown_ref, &child_vis, 10);

    // deleted=2 (child 0 gap + child 1 variable indel), present=0
    // deleted == n_children: resolved as consensus gap, removed from variable_indel
    assert!(result.variable_indel.is_empty());
    assert_eq!(result.resolved_gaps, vec![(2, 4)]);
  }

  #[test]
  fn test_fitch_indel_forward_variable_parent_has_sequence() {
    // Variable indel (majority deleted), parent has sequence: parent wins → no deletion, no gap.
    let mut variable_indel = BTreeMap::new();
    variable_indel.insert((2, 5), Deletion { deleted: 2, present: 1 });

    let parent_gaps: Vec<(usize, usize)> = vec![];
    let parent_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let node_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let node_gaps: Vec<(usize, usize)> = vec![];

    let (indels, node_gaps) = resolve_indels_forward(&variable_indel, &node_gaps, &[], &parent_gaps, &parent_seq, &node_seq);

    assert!(indels.is_empty());
    assert!(node_gaps.is_empty());
  }

  #[test]
  fn test_fitch_indel_forward_variable_parent_has_gap() {
    // Variable indel (majority present), parent has gap: parent wins → node inherits gap, no insertion.
    let mut variable_indel = BTreeMap::new();
    variable_indel.insert((2, 5), Deletion { deleted: 1, present: 2 });

    let parent_gaps = vec![(2, 5)];
    let parent_seq = Seq::try_from_str("AC---CGTAC").unwrap();
    let node_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let node_gaps: Vec<(usize, usize)> = vec![];

    let (indels, node_gaps) = resolve_indels_forward(&variable_indel, &node_gaps, &[], &parent_gaps, &parent_seq, &node_seq);

    assert!(indels.is_empty());
    assert_eq!(node_gaps, vec![(2, 5)]);
  }

  #[test]
  fn test_fitch_indel_forward_consensus_deletion() {
    let variable_indel = BTreeMap::new();
    let parent_gaps: Vec<(usize, usize)> = vec![];
    let parent_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let node_seq = Seq::try_from_str("AC---CGTAC").unwrap();
    let node_gaps = vec![(2, 5)];

    let (indels, node_gaps) = resolve_indels_forward(&variable_indel, &node_gaps, &[], &parent_gaps, &parent_seq, &node_seq);

    assert_eq!(indels.len(), 1);
    assert!(indels[0].deletion);
    assert_eq!(indels[0].range, (2, 5));
    assert_eq!(indels[0].seq, Seq::try_from_str("GTA").unwrap());
  }

  #[test]
  fn test_fitch_indel_forward_consensus_insertion() {
    let variable_indel = BTreeMap::new();
    let parent_gaps = vec![(2, 5)];
    let parent_seq = Seq::try_from_str("AC---CGTAC").unwrap();
    let node_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let node_gaps: Vec<(usize, usize)> = vec![];

    let (indels, node_gaps) = resolve_indels_forward(&variable_indel, &node_gaps, &[], &parent_gaps, &parent_seq, &node_seq);

    assert_eq!(indels.len(), 1);
    assert!(!indels[0].deletion);
    assert_eq!(indels[0].range, (2, 5));
    assert_eq!(indels[0].seq, Seq::try_from_str("GTA").unwrap());
  }

  #[test]
  fn test_fitch_indel_forward_no_indels() {
    let variable_indel = BTreeMap::new();
    let parent_gaps: Vec<(usize, usize)> = vec![];
    let parent_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let node_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let node_gaps: Vec<(usize, usize)> = vec![];

    let (indels, node_gaps) = resolve_indels_forward(&variable_indel, &node_gaps, &[], &parent_gaps, &parent_seq, &node_seq);

    assert!(indels.is_empty());
  }
}
