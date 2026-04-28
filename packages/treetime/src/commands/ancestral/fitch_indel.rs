use crate::representation::payload::sparse::Deletion;
use crate::seq::indel::InDel;
use std::collections::BTreeMap;
use treetime_primitives::Seq;
use treetime_utils::interval::range_complement::range_complement;
use treetime_utils::interval::range_difference::range_difference;
use treetime_utils::interval::range_intersection::range_intersection;

/// Resolve indels during the backward pass: identify positions where children
/// disagree on gap presence.
///
/// Operates on gap range lists and integer counters. Independent of sparse/dense
/// types - usable by both partition kinds.
///
/// Three steps, matching the sparse Fitch backward logic:
/// 1. For each child, intersect child gaps with complement of consensus gaps
///    (intersection of all children's gaps) to find disagreement ranges.
/// 2. Propagate variable indels from children.
/// 3. Caller is responsible for collapsing ranges where all children agree on gap.
///
/// Returns the variable_indel map. The caller must filter out entries where
/// `deleted == n_children` and push those ranges back to the consensus gaps.
pub fn resolve_indels_backward(
  child_gaps: &[Vec<(usize, usize)>],
  child_variable_indels: &[&BTreeMap<(usize, usize), Deletion>],
  length: usize,
) -> BTreeMap<(usize, usize), Deletion> {
  let n_children = child_gaps.len();
  let consensus_gaps = range_intersection(child_gaps);
  let non_gap = range_complement(&[(0, length)], &[consensus_gaps]);

  let mut variable_indel = BTreeMap::new();

  // Step 1: find gap ranges in children that are not in consensus (disagreement)
  for child_gap in child_gaps {
    for r in range_intersection(&[non_gap.clone(), child_gap.clone()]) {
      let indel = variable_indel.entry(r).or_insert_with(|| Deletion {
        deleted: 0,
        present: n_children,
      });
      indel.deleted += 1;
      indel.present -= 1;
    }
  }

  // Step 2: propagate variable indels from children
  for child_vi in child_variable_indels {
    for r in child_vi.keys() {
      if let Some(indel) = variable_indel.get_mut(r) {
        indel.deleted += 1;
        indel.present -= 1;
      }
    }
  }

  variable_indel
}

/// Resolve variable indels during the forward pass using parent context.
///
/// For each variable indel range, uses parent gap state to break ties.
/// Also detects consensus gap differences between node and parent (non-variable
/// deletions/insertions).
///
/// Returns the list of resolved InDel entries for the edge.
///
/// Dense partitions do not need `InDel::seq` content for inserted sequences
/// because the full sequence is available from the dense reconstruction.
/// We populate it for consistency with sparse indels.
pub fn resolve_indels_forward(
  variable_indel: &BTreeMap<(usize, usize), Deletion>,
  node_gaps: &mut Vec<(usize, usize)>,
  parent_gaps: &[(usize, usize)],
  parent_sequence: &Seq,
  node_sequence: &Seq,
) -> Vec<InDel> {
  let mut indels = Vec::new();

  // Process variable indels using parent context for tiebreaking
  for (r, indel) in variable_indel {
    let gap_in_parent = if parent_gaps.contains(r) { 1 } else { 0 };
    if indel.deleted + gap_in_parent > indel.present {
      node_gaps.push(*r);
      if gap_in_parent == 0 {
        let indel = InDel::del(*r, &parent_sequence[r.0..r.1]);
        indels.push(indel);
      }
    } else if gap_in_parent > 0 {
      let indel = InDel::ins(*r, &node_sequence[r.0..r.1]);
      indels.push(indel);
    }
  }

  // Process consensus gaps in node that are not in parent (deletions)
  for r in range_difference(node_gaps, parent_gaps) {
    if variable_indel.contains_key(&r) {
      continue;
    }
    let indel = InDel::del(r, &node_sequence[r.0..r.1]);
    indels.push(indel);
  }

  // Process gaps in parent that are not in node (insertions)
  for r in range_difference(parent_gaps, node_gaps) {
    if variable_indel.contains_key(&r) {
      continue;
    }
    let indel = InDel::ins(r, &node_sequence[r.0..r.1]);
    indels.push(indel);
  }

  indels
}
