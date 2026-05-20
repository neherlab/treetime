use crate::representation::payload::sparse::Deletion;
use crate::seq::indel::InDel;
use std::collections::BTreeMap;
use itertools::Itertools;
use treetime_primitives::Seq;
use treetime_utils::interval::range_difference::range_difference;
use treetime_utils::interval::range_intersection::{range_intersection, range_intersection_iter};
use treetime_utils::interval::range_union::range_union_iter;

/// Collect sorted, deduplicated breakpoints from an iterator of `(lo, hi)` ranges.
fn collect_breakpoints(ranges: impl Iterator<Item = (usize, usize)>) -> Vec<usize> {
  let mut bp: Vec<usize> = ranges.flat_map(|(lo, hi)| [lo, hi]).collect();
  bp.sort_unstable();
  bp.dedup();
  bp
}

/// Returns true when `[lo, hi)` lies entirely within some range of `ranges`.
fn interval_in(ranges: &[(usize, usize)], lo: usize, hi: usize) -> bool {
  ranges.iter().any(|&(a, b)| a <= lo && hi <= b)
}

/// Returns true when `[lo, hi)` lies entirely within some key of `vi`.
fn interval_in_vi(vi: &BTreeMap<(usize, usize), Deletion>, lo: usize, hi: usize) -> bool {
  vi.keys().any(|&(a, b)| a <= lo && hi <= b)
}


/// Ranges describing a parent node's gap/unknown/non_char state, derived from its children.
pub struct NodeRanges {
  /// Positions where all children lack data (gap or unknown).
  pub non_char: Vec<(usize, usize)>,
  /// Positions in non_char where no child has a gap.
  pub unknown: Vec<(usize, usize)>,
}

/// Compute the parent node's non_char, consensus_gaps, and unknown ranges from child data.
pub fn compute_node_ranges(
  child_non_chars: &[&Vec<(usize, usize)>],
  child_gaps: &[&Vec<(usize, usize)>],
) -> NodeRanges {
  let non_char = range_intersection_iter(child_non_chars.iter().copied()).collect_vec();
  let gap_union = range_union_iter(child_gaps.iter().copied()).collect_vec();
  let consensus_gaps = range_intersection(&[non_char.clone(), gap_union]);
  let unknown = range_difference(&non_char, &consensus_gaps);
  NodeRanges { non_char, unknown }
}

pub struct IndelsBackward {
  pub variable_indel: BTreeMap<(usize, usize), Deletion>,
  pub resolved_gaps: Vec<(usize, usize)>,
}

/// Resolve indels during the backward pass: identify positions where children
/// disagree on gap presence.
///
/// Returns variable indels (unresolved disagreements) and resolved gaps (ranges
/// where all children agree on gap).
pub fn resolve_indels_backward(
  child_gaps: &[&Vec<(usize, usize)>],
  child_unknown: &[&Vec<(usize, usize)>],
  child_variable_indels: &[&BTreeMap<(usize, usize), Deletion>],
  length: usize,
) -> IndelsBackward {
  let n_children = child_gaps.len();

  if n_children == 1 {
    return IndelsBackward {
      variable_indel: child_variable_indels[0].clone(),
      resolved_gaps: child_gaps[0].clone(),
    };
  }

  // Collect all interval breakpoints from every child source so we can sweep
  // over non-overlapping atomic sub-intervals, one entry per interval.
  let breakpoints = collect_breakpoints(
    std::iter::once((0, length))
      .chain(child_gaps.iter().flat_map(|g| g.iter().copied()))
      .chain(child_unknown.iter().flat_map(|u| u.iter().copied()))
      .chain(child_variable_indels.iter().flat_map(|vi| vi.keys().copied())),
  );

  // Sweep over atomic intervals, classifying each by child gap evidence.
  // Collect and merge intervals as resolved consensus gaps, unresolved variable indels, or skip for no deletion
  let mut variable_indel: BTreeMap<(usize, usize), Deletion> = BTreeMap::new();
  let mut resolved_gaps: Vec<(usize, usize)> = Vec::new();

  for w in breakpoints.windows(2) {
    let (lo, hi) = (w[0], w[1]);

    // Classify each child's state at [lo, hi).
    // Priority: gapped > variable_indel > unknown > char.
    let mut n_gapped: usize = 0;
    let mut n_variable: usize = 0;
    let mut n_unknown: usize = 0;
    for i in 0..n_children {
      if interval_in(&child_gaps[i], lo, hi) {
        n_gapped += 1;
      } else if interval_in_vi(child_variable_indels[i], lo, hi) {
        n_variable += 1;
      } else if interval_in(&child_unknown[i], lo, hi) {
        n_unknown += 1;
      }
    }

    // No gap evidence at this interval, but at least one child as sequence --> no gap at parent.
    // OR All children unknown --> no-gap evidence.
    if (n_gapped == 0 && (n_variable + n_unknown < n_children)) || n_unknown == n_children {
      continue;
    }

    let n_no_seq = n_gapped + n_unknown;
    if n_gapped > 0 && (n_no_seq + n_variable == n_children) { // Evidence for gap and all children compatible with gap --> resolved gap
      if let Some(last) = resolved_gaps.last_mut() {
        if last.1 == lo {
          last.1 = hi;
        } else {
          resolved_gaps.push((lo, hi));
        }
      } else {
        resolved_gaps.push((lo, hi));
      }
    } else { // if n_gapped==0, then n_variable + n_unknown == n_children and n_variable>0 --> can't resolve, keep variable.
             // If n_gapped>0, then n_no_seq + n_variable < n_children --> there are children with and without sequence --> variable
      let del = Deletion { deleted: n_no_seq, present: n_children - n_no_seq };
      variable_indel.insert((lo, hi), del);
    }
  }

  IndelsBackward {
    variable_indel,
    resolved_gaps,
  }
}

/// Resolve variable indels during the forward pass using parent context.
///
/// For each variable indel range, uses parent gap state to break ties.
/// Also detects consensus gap differences between node and parent (non-variable
/// deletions/insertions).
///
/// Returns the list of resolved InDel entries for the edge and the final gap set of the node.
///
/// Dense partitions do not need `InDel::seq` content for inserted sequences
/// because the full sequence is available from the dense reconstruction.
/// We populate it for consistency with sparse indels.
pub fn resolve_indels_forward(
  variable_indel: &BTreeMap<(usize, usize), Deletion>,
  node_gaps: &[(usize, usize)],
  node_non_char: &[(usize, usize)],
  parent_gaps: &[(usize, usize)],
  parent_sequence: &Seq,
  node_sequence: &Seq,
) -> (Vec<InDel>, Vec<(usize, usize)>) {
  // Collect all breakpoints from every relevant source.
  let breakpoints = collect_breakpoints(
    parent_gaps.iter().copied()
      .chain(node_gaps.iter().copied())
      .chain(node_non_char.iter().copied())
      .chain(variable_indel.keys().copied()),
  );

  let mut new_node_gaps: Vec<(usize, usize)> = Vec::new();
  let mut deletions: Vec<InDel> = Vec::new();
  let mut insertions: Vec<InDel> = Vec::new();

  for w in breakpoints.windows(2) {
    // the semantics here refer to an interval [lo, hi) and generally to the presence of gap
    let (lo, hi) = (w[0], w[1]);
    let in_parent = interval_in(parent_gaps, lo, hi);
    let in_node = interval_in(node_gaps, lo, hi);
    let in_non_char = interval_in(node_non_char, lo, hi);
    let in_vi = interval_in_vi(variable_indel, lo, hi);

    if in_node && !in_parent {
      // Node has gap, parent doesn't: deletion on this edge.
      if let Some(last) = deletions.last_mut().filter(|d| d.range.1 == lo) {
        last.range.1 = hi;
        last.seq.extend(parent_sequence[lo..hi].iter().copied());
      } else {
        deletions.push(InDel::del((lo, hi), &parent_sequence[lo..hi]));
      }
      // Gap persists in the node.
      if let Some(last) = new_node_gaps.last_mut().filter(|l| l.1 == lo) {
        last.1 = hi;
      } else {
        new_node_gaps.push((lo, hi));
      }
    } else if in_parent && !in_node && !in_non_char && !in_vi {
      // Parent has gap, node is not gap-compatible: insertion on this edge.
      if let Some(last) = insertions.last_mut().filter(|i| i.range.1 == lo) {
        last.range.1 = hi;
        last.seq.extend(node_sequence[lo..hi].iter().copied());
      } else {
        insertions.push(InDel::ins((lo, hi), &node_sequence[lo..hi]));
      }
    } else if in_parent {
      // Parent has gap; node is gapped, non_char, or variable → inherit parent gap.
      if let Some(last) = new_node_gaps.last_mut().filter(|l| l.1 == lo) {
        last.1 = hi;
      } else {
        new_node_gaps.push((lo, hi));
      }
    }
    // !in_parent && !in_node: no gap on either side, nothing to do.
  }

  (deletions.into_iter().chain(insertions).collect(), new_node_gaps)
}
