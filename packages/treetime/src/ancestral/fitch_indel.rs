use crate::representation::payload::sparse::Deletion;
use crate::seq::indel::InDel;
use std::collections::BTreeMap;
use itertools::Itertools;
use log::debug;
use treetime_primitives::Seq;
use treetime_utils::interval::range_difference::range_difference;
use treetime_utils::interval::range_intersection::{range_intersection, range_intersection_iter};
use treetime_utils::interval::range_union::range_union_iter;

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
  /// Positions in non_char where at least one child has a gap (gap beats unknown).
  pub consensus_gaps: Vec<(usize, usize)>,
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
  NodeRanges { non_char, consensus_gaps, unknown }
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
  consensus_gaps: Vec<(usize, usize)>,
  child_gaps: &[&Vec<(usize, usize)>],
  child_unknown: &[&Vec<(usize, usize)>],
  child_variable_indels: &[&BTreeMap<(usize, usize), Deletion>],
  length: usize,
) -> IndelsBackward {
  debug!(
    "resolve_indels_backward: length={length}, consensus_gaps={consensus_gaps:?}, child_gaps={child_gaps:?}, child_unknown={child_unknown:?}, child_variable_indels={child_variable_indels:?}"
  );
  let n_children = child_gaps.len();

  // Collect all interval breakpoints from every child source so we can sweep
  // over non-overlapping atomic sub-intervals.
  let mut breakpoints: Vec<usize> = vec![0, length];
  for &(lo, hi) in &consensus_gaps {
    breakpoints.push(lo);
    breakpoints.push(hi);
  }
  for child_gap in child_gaps {
    for &(lo, hi) in *child_gap {
      breakpoints.push(lo);
      breakpoints.push(hi);
    }
  }
  for child_unk in child_unknown {
    for &(lo, hi) in *child_unk {
      breakpoints.push(lo);
      breakpoints.push(hi);
    }
  }
  for child_vi in child_variable_indels {
    for &(lo, hi) in child_vi.keys() {
      breakpoints.push(lo);
      breakpoints.push(hi);
    }
  }
  breakpoints.sort_unstable();
  breakpoints.dedup();

  // Sweep over atomic intervals, accumulating deletion counts for each.
  // We merge consecutive intervals that carry identical counts as we go.
  let mut variable_indel: BTreeMap<(usize, usize), Deletion> = BTreeMap::new();
  let mut resolved_gaps: Vec<(usize, usize)> = Vec::new();

  // `pending` holds the current merged interval being built.
  let mut pending: Option<((usize, usize), Deletion)> = None;

  let flush = |pending: Option<((usize, usize), Deletion)>,
                   variable_indel: &mut BTreeMap<(usize, usize), Deletion>,
                   resolved_gaps: &mut Vec<(usize, usize)>,
                   n_children: usize| {
    if let Some(((lo, hi), del)) = pending {
      if del.deleted == n_children {
        resolved_gaps.push((lo, hi));
      } else {
        variable_indel.insert((lo, hi), del);
      }
    }
  };

  for w in breakpoints.windows(2) {
    let (lo, hi) = (w[0], w[1]);

    // Skip intervals that fall inside the consensus gaps.
    if interval_in(&consensus_gaps, lo, hi) {
      flush(pending.take(), &mut variable_indel, &mut resolved_gaps, n_children);
      continue;
    }

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

    // No gap evidence at this interval: nothing to record.
    if n_gapped == 0 && n_variable == 0 {
      flush(pending.take(), &mut variable_indel, &mut resolved_gaps, n_children);
      continue;
    }

    // Unknowns and child variable indels lean toward "deleted" (gap-beats-unknown rule).
    let n_deleted = n_gapped + n_unknown + n_variable;
    let del = Deletion { deleted: n_deleted, present: n_children - n_deleted };

    // Merge with the pending interval if it is adjacent and has the same counts.
    match &mut pending {
      Some(((_, cur_hi), cur_del)) if *cur_hi == lo && cur_del.deleted == del.deleted => {
        *cur_hi = hi;
      },
      _ => {
        flush(pending.take(), &mut variable_indel, &mut resolved_gaps, n_children);
        pending = Some(((lo, hi), del));
      },
    }
  }
  flush(pending, &mut variable_indel, &mut resolved_gaps, n_children);

  debug!(
    "resolve_indels_backward result: variable_indel={variable_indel:?}, resolved_gaps={resolved_gaps:?}"
  );
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
/// Returns the list of resolved InDel entries for the edge.
///
/// Dense partitions do not need `InDel::seq` content for inserted sequences
/// because the full sequence is available from the dense reconstruction.
/// We populate it for consistency with sparse indels.
pub fn resolve_indels_forward(
  variable_indel: &BTreeMap<(usize, usize), Deletion>,
  node_gaps: &mut Vec<(usize, usize)>,
  node_non_char: &[(usize, usize)],
  parent_gaps: &[(usize, usize)],
  parent_sequence: &Seq,
  node_sequence: &Seq,
) -> Vec<InDel> {
  let mut indels = Vec::new();
  debug!("resolve_indels_forward");
  // Process variable indels using parent context for tiebreaking
  for (r, indel) in variable_indel {
    let gap_in_parent = if parent_gaps.contains(r) { 1 } else { 0 };
    debug!(
      "resolve_indels_forward variable_indel: range={r:?}, indel={indel:?}, gap_in_parent={gap_in_parent}, parent_seq={:?}, node_seq={:?}",
      Seq::from(&parent_sequence[r.0..r.1]).as_str(),
      Seq::from(&node_sequence[r.0..r.1]).as_str(),
    );
    if indel.deleted + gap_in_parent > indel.present {
      node_gaps.push(*r);
      if gap_in_parent == 0 {
        let indel = InDel::del(*r, &parent_sequence[r.0..r.1]);
        indels.push(indel);
      }
    } else if gap_in_parent > 0 && !interval_in(node_non_char, r.0, r.1) {
      let indel = InDel::ins(*r, &node_sequence[r.0..r.1]);
      indels.push(indel);
    }
  }

  // Process consensus gaps in node that are not in parent (deletions)
  for r in range_difference(node_gaps, parent_gaps) {
    if variable_indel.contains_key(&r) {
      continue;
    }
    let indel = InDel::del(r, &parent_sequence[r.0..r.1]);
    indels.push(indel);
  }

  // Process gaps in parent that are not in node (insertions).
  // Skip ranges where the node has no known sequence (non_char): those positions
  // were filled from the parent during the forward pass and do not represent
  // an actual insertion of new sequence.
  for r in range_difference(parent_gaps, node_gaps) {
    if variable_indel.contains_key(&r) || interval_in(node_non_char, r.0, r.1) {
      continue;
    }
    let indel = InDel::ins(r, &node_sequence[r.0..r.1]);
    indels.push(indel);
  }

  indels
}
