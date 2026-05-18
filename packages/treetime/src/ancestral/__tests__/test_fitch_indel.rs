#[cfg(test)]
mod tests {
  use crate::ancestral::fitch_indel::{resolve_indels_backward, resolve_indels_forward};
  use crate::representation::payload::sparse::Deletion;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_primitives::Seq;

  #[test]
  fn test_fitch_indel_backward_no_disagreement() {
    let child_gaps = vec![vec![(2, 4)], vec![(2, 4)]];
    let consensus_gaps = vec![(2, 4)];
    let child_unknown = vec![vec![], vec![]];
    let empty0 = BTreeMap::new();
    let empty1 = BTreeMap::new();
    let child_vis: Vec<&BTreeMap<(usize, usize), Deletion>> = vec![&empty0, &empty1];
    let result = resolve_indels_backward(consensus_gaps, &child_gaps, &child_unknown, &child_vis, 10).variable_indel;
    assert!(result.is_empty(), "No disagreement when all children have same gaps");
  }

  #[test]
  fn test_fitch_indel_backward_one_child_has_gap() {
    let child_gaps = vec![vec![(2, 4)], vec![]];
    let consensus_gaps = vec![];
    let child_unknown = vec![vec![], vec![]];
    let empty0 = BTreeMap::new();
    let empty1 = BTreeMap::new();
    let child_vis: Vec<&BTreeMap<(usize, usize), Deletion>> = vec![&empty0, &empty1];
    let result = resolve_indels_backward(consensus_gaps, &child_gaps, &child_unknown, &child_vis, 10).variable_indel;
    assert_eq!(result.len(), 1);
    let del = &result[&(2, 4)];
    assert_eq!(del.deleted, 1);
    assert_eq!(del.present, 1);
  }

  #[test]
  fn test_fitch_indel_backward_all_children_deleted() {
    let child_gaps = vec![vec![(2, 4)], vec![(2, 4)], vec![(2, 4)]];
    let consensus_gaps = vec![(2, 4)];
    let child_unknown = vec![vec![], vec![], vec![]];
    let empty0 = BTreeMap::new();
    let empty1 = BTreeMap::new();
    let empty2 = BTreeMap::new();
    let child_vis: Vec<&BTreeMap<(usize, usize), Deletion>> = vec![&empty0, &empty1, &empty2];

    // consensus gaps = intersection = [(2,4)], so this range is skipped.
    // No disagreement ranges produced.
    let result = resolve_indels_backward(consensus_gaps, &child_gaps, &child_unknown, &child_vis, 10).variable_indel;
    assert!(result.is_empty(), "All children agree on gap, no variable indels");
  }

  #[test]
  fn test_fitch_indel_backward_propagates_child_variable_indels() {
    // Child 0 has gap at (2,4), child 1 does not.
    // Child 1 has variable indel at (2,4).
    let child_gaps = vec![vec![(2, 4)], vec![]];
    let consensus_gaps = vec![];
    let child_unknown = vec![vec![], vec![]];
    let empty0 = BTreeMap::new();
    let mut child1_vi = BTreeMap::new();
    child1_vi.insert((2, 4), Deletion { deleted: 1, present: 1 });
    let child_vis: Vec<&BTreeMap<(usize, usize), Deletion>> = vec![&empty0, &child1_vi];
    let result = resolve_indels_backward(consensus_gaps, &child_gaps, &child_unknown, &child_vis, 10);

    // deleted=2 (child 0 gap + child 1 variable indel), present=0
    // deleted == n_children: resolved as consensus gap, removed from variable_indel
    assert!(result.variable_indel.is_empty());
    assert_eq!(result.resolved_gaps, vec![(2, 4)]);
  }

  #[test]
  fn test_fitch_indel_forward_deletion_from_variable() {
    let mut variable_indel = BTreeMap::new();
    variable_indel.insert((2, 5), Deletion { deleted: 2, present: 1 });

    let parent_gaps: Vec<(usize, usize)> = vec![];
    let parent_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let node_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let mut node_gaps: Vec<(usize, usize)> = vec![];

    let indels = resolve_indels_forward(&variable_indel, &mut node_gaps, &[], &parent_gaps, &parent_seq, &node_seq);

    assert_eq!(indels.len(), 1);
    assert!(indels[0].deletion);
    assert_eq!(indels[0].range, (2, 5));
    assert_eq!(indels[0].seq, Seq::try_from_str("GTA").unwrap());
    assert_eq!(node_gaps, vec![(2, 5)]);
  }

  #[test]
  fn test_fitch_indel_forward_insertion_from_variable() {
    let mut variable_indel = BTreeMap::new();
    variable_indel.insert((2, 5), Deletion { deleted: 1, present: 2 });

    let parent_gaps = vec![(2, 5)];
    let parent_seq = Seq::try_from_str("AC---CGTAC").unwrap();
    let node_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let mut node_gaps: Vec<(usize, usize)> = vec![];

    let indels = resolve_indels_forward(&variable_indel, &mut node_gaps, &[], &parent_gaps, &parent_seq, &node_seq);

    assert_eq!(indels.len(), 1);
    assert!(!indels[0].deletion);
    assert_eq!(indels[0].range, (2, 5));
    assert_eq!(indels[0].seq, Seq::try_from_str("GTA").unwrap());
  }

  #[test]
  fn test_fitch_indel_forward_consensus_deletion() {
    let variable_indel = BTreeMap::new();
    let parent_gaps: Vec<(usize, usize)> = vec![];
    let parent_seq = Seq::try_from_str("ACGTACGTAC").unwrap();
    let node_seq = Seq::try_from_str("AC---CGTAC").unwrap();
    let mut node_gaps = vec![(2, 5)];

    let indels = resolve_indels_forward(&variable_indel, &mut node_gaps, &[], &parent_gaps, &parent_seq, &node_seq);

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
    let mut node_gaps: Vec<(usize, usize)> = vec![];

    let indels = resolve_indels_forward(&variable_indel, &mut node_gaps, &[], &parent_gaps, &parent_seq, &node_seq);

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
    let mut node_gaps: Vec<(usize, usize)> = vec![];

    let indels = resolve_indels_forward(&variable_indel, &mut node_gaps, &[], &parent_gaps, &parent_seq, &node_seq);

    assert!(indels.is_empty());
  }
}
