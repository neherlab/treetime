use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use std::fmt;
use treetime_primitives::Seq;
use treetime_utils::interval::range_difference::range_difference;
use treetime_utils::interval::range_intersection::{range_intersection, range_intersection_iter};
use treetime_utils::interval::range_union::range_union_iter;

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct InDel {
  pub range: (usize, usize),
  pub seq: Seq,
  pub deletion: bool, // deletion if True, insertion if False
}

impl InDel {
  pub fn del(range: (usize, usize), seq: impl Into<Seq>) -> Self {
    Self::new(range, seq, true)
  }

  pub fn ins(range: (usize, usize), seq: impl Into<Seq>) -> Self {
    Self::new(range, seq, false)
  }

  pub fn new(range: (usize, usize), seq: impl Into<Seq>, deletion: bool) -> Self {
    let seq = seq.into();
    assert!(range.0 <= range.1);
    assert_eq!(seq.len(), range.1 - range.0);
    Self { range, seq, deletion }
  }

  /// Invert the indel direction by toggling deletion flag
  pub fn invert(&mut self) {
    self.deletion = !self.deletion;
  }
}

impl fmt::Display for InDel {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    let delta_str = if self.deletion {
      format!("{} -> {}", &self.seq.as_str(), "-".repeat(self.seq.len()))
    } else {
      format!("{} -> {}", "-".repeat(self.seq.len()), &self.seq.as_str())
    };
    write!(f, "{}--{}: {}", self.range.0, self.range.1, delta_str)
  }
}

/// Compose indels from two consecutive edges into net indels.
///
/// When collapsing node B between parent A and child C, indels on edges A->B
/// (`parent_indels`) and B->C (`child_indels`) are composed to produce the net
/// A->C indels. This parallels `compose_substitutions()` for substitutions.
///
/// Composition rules by case:
/// - Non-overlapping: both pass through unchanged
/// - Adjacent deletions: merged into a single deletion spanning both ranges
/// - Overlapping deletions (same or partial): unified using parent's seq content
/// - Overlapping insertions (same range): child's seq wins
/// - Deletion then insertion (same range, same content): cancel
/// - Deletion then insertion (same range, different content): both kept
///
/// Both input slices must be sorted by `range.0` (ties broken by `range.1`).
/// The output is sorted by `range.0`.
pub fn compose_indels(parent_indels: &[InDel], child_indels: &[InDel]) -> Vec<InDel> {
  debug_assert!(
    parent_indels
      .windows(2)
      .all(|w| matches!(w, [a, b] if a.range <= b.range)),
    "parent_indels not sorted by range"
  );
  debug_assert!(
    child_indels
      .windows(2)
      .all(|w| matches!(w, [a, b] if a.range <= b.range)),
    "child_indels not sorted by range"
  );

  if parent_indels.is_empty() {
    return child_indels.to_vec();
  }
  if child_indels.is_empty() {
    return parent_indels.to_vec();
  }

  let mut result = Vec::with_capacity(parent_indels.len() + child_indels.len());
  let mut pi = 0;
  let mut ci = 0;

  while pi < parent_indels.len() && ci < child_indels.len() {
    let p = &parent_indels[pi];
    let c = &child_indels[ci];

    let overlap_start = p.range.0.max(c.range.0);
    let overlap_end = p.range.1.min(c.range.1);
    let overlaps = overlap_start < overlap_end;
    let adjacent = !overlaps && p.range.1 == c.range.0 && p.deletion && c.deletion;

    if !overlaps && !adjacent {
      // Case 7: no interaction
      if p.range.0 < c.range.0 || (p.range.0 == c.range.0 && p.range.1 <= c.range.1) {
        result.push(p.clone());
        pi += 1;
      } else {
        result.push(c.clone());
        ci += 1;
      }
    } else if adjacent {
      // Case 2: adjacent deletions - merge into one
      let merged_seq: Seq = p
        .seq
        .as_slice()
        .iter()
        .chain(c.seq.as_slice().iter())
        .copied()
        .collect();
      result.push(InDel::del((p.range.0, c.range.1), merged_seq));
      pi += 1;
      ci += 1;
    } else {
      // Overlapping cases - both pointers always advance
      match (p.deletion, c.deletion) {
        (true, true) => {
          // Cases 1, 6: overlapping deletions - use parent's seq for overlap region
          let merged_start = p.range.0.min(c.range.0);
          let merged_end = p.range.1.max(c.range.1);

          let mut merged_seq = Seq::with_capacity(merged_end - merged_start);

          if p.range.0 < c.range.0 {
            let prefix_len = c.range.0 - p.range.0;
            merged_seq.extend(p.seq.as_slice()[..prefix_len].iter().copied());
          } else if c.range.0 < p.range.0 {
            let prefix_len = p.range.0 - c.range.0;
            merged_seq.extend(c.seq.as_slice()[..prefix_len].iter().copied());
          }

          // Overlap region: parent's seq has the original content
          let p_overlap_offset = overlap_start - p.range.0;
          let p_overlap_len = overlap_end - overlap_start;
          merged_seq.extend(
            p.seq.as_slice()[p_overlap_offset..p_overlap_offset + p_overlap_len]
              .iter()
              .copied(),
          );

          if p.range.1 > c.range.1 {
            let suffix_offset = overlap_end - p.range.0;
            merged_seq.extend(p.seq.as_slice()[suffix_offset..].iter().copied());
          } else if c.range.1 > p.range.1 {
            let suffix_offset = overlap_end - c.range.0;
            merged_seq.extend(c.seq.as_slice()[suffix_offset..].iter().copied());
          }

          result.push(InDel::del((merged_start, merged_end), merged_seq));
        },
        (false, false) => {
          // Case 3: overlapping insertions - child's seq wins
          result.push(c.clone());
        },
        (true, false) | (false, true) => {
          if !(p.range == c.range && p.seq == c.seq) {
            // Emit in sorted order by range.0 to maintain output sortedness
            if p.range.0 <= c.range.0 {
              result.push(p.clone());
              result.push(c.clone());
            } else {
              result.push(c.clone());
              result.push(p.clone());
            }
          }
        },
      }
      pi += 1;
      ci += 1;
    }
  }

  result.extend_from_slice(&parent_indels[pi..]);
  result.extend_from_slice(&child_indels[ci..]);

  // Merge adjacent deletions in a final pass (handles chains from multi-step composition)
  merge_adjacent_deletions(result)
}

fn merge_adjacent_deletions(indels: Vec<InDel>) -> Vec<InDel> {
  if indels.len() <= 1 {
    return indels;
  }

  let mut merged: Vec<InDel> = Vec::with_capacity(indels.len());
  for indel in indels {
    #[allow(clippy::suspicious_operation_groupings)]
    let should_merge = merged
      .last()
      .is_some_and(|prev: &InDel| prev.deletion && indel.deletion && prev.range.1 == indel.range.0);
    if should_merge {
      let prev = merged.last_mut().unwrap();
      prev.range.1 = indel.range.1;
      prev.seq.extend(indel.seq.as_slice().iter().copied());
    } else {
      merged.push(indel);
    }
  }
  merged
}

/// Sort indels by range for use as input to `compose_indels`.
pub fn sort_indels(indels: &mut [InDel]) {
  indels.sort_by_key(|i| i.range);
}

pub struct IndelsBackward {
  pub variable_indel: BTreeSet<(usize, usize)>,
  pub resolved_gaps: Vec<(usize, usize)>,
}

pub struct NodeRanges {
  pub non_char: Vec<(usize, usize)>,
  pub unknown: Vec<(usize, usize)>,
}

fn collect_breakpoints(ranges: impl Iterator<Item = (usize, usize)>) -> Vec<usize> {
  let mut bp: Vec<usize> = ranges.flat_map(|(lo, hi)| <[usize; 2]>::from((lo, hi))).collect();
  bp.sort_unstable();
  bp.dedup();
  bp
}

fn interval_in(ranges: &[(usize, usize)], lo: usize, hi: usize) -> bool {
  ranges.iter().any(|&(a, b)| a <= lo && hi <= b)
}

fn interval_in_vi(vi: &BTreeSet<(usize, usize)>, lo: usize, hi: usize) -> bool {
  vi.iter().any(|&(a, b)| a <= lo && hi <= b)
}

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

pub fn resolve_indels_backward(
  child_gaps: &[&Vec<(usize, usize)>],
  child_unknown: &[&Vec<(usize, usize)>],
  child_variable_indels: &[&BTreeSet<(usize, usize)>],
  length: usize,
) -> IndelsBackward {
  let n_children = child_gaps.len();

  if n_children == 1 {
    return IndelsBackward {
      variable_indel: child_variable_indels[0].clone(),
      resolved_gaps: child_gaps[0].clone(),
    };
  }

  let breakpoints = collect_breakpoints(
    std::iter::once((0, length))
      .chain(child_gaps.iter().flat_map(|g| g.iter().copied()))
      .chain(child_unknown.iter().flat_map(|u| u.iter().copied()))
      .chain(child_variable_indels.iter().flat_map(|vi| vi.iter().copied())),
  );

  let mut variable_indel: BTreeSet<(usize, usize)> = BTreeSet::new();
  let mut resolved_gaps: Vec<(usize, usize)> = Vec::new();

  for (lo, hi) in breakpoints.iter().copied().tuple_windows() {
    // Classify each child's state at [lo, hi).
    // Priority: gapped > variable_indel > unknown > char.
    let mut n_gapped: usize = 0;
    let mut n_variable: usize = 0;
    let mut n_unknown: usize = 0;
    for i in 0..n_children {
      if interval_in(child_gaps[i], lo, hi) {
        n_gapped += 1;
      } else if interval_in_vi(child_variable_indels[i], lo, hi) {
        n_variable += 1;
      } else if interval_in(child_unknown[i], lo, hi) {
        n_unknown += 1;
      }
    }

    // No gap evidence, or all children unknown: skip.
    if (n_gapped == 0 && (n_variable + n_unknown < n_children)) || n_unknown == n_children {
      continue;
    }

    let n_no_seq = n_gapped + n_unknown;
    if n_gapped > 0 && (n_no_seq + n_variable == n_children) {
      // All children compatible with gap and at least one has explicit gap: resolved.
      if let Some(last) = resolved_gaps.last_mut() {
        if last.1 == lo {
          last.1 = hi;
        } else {
          resolved_gaps.push((lo, hi));
        }
      } else {
        resolved_gaps.push((lo, hi));
      }
    } else {
      variable_indel.insert((lo, hi));
    }
  }

  IndelsBackward {
    variable_indel,
    resolved_gaps,
  }
}

pub fn resolve_indels_forward(
  variable_indel: &BTreeSet<(usize, usize)>,
  node_gaps: &[(usize, usize)],
  node_non_char: &[(usize, usize)],
  parent_gaps: &[(usize, usize)],
  parent_sequence: &Seq,
  node_sequence: &Seq,
) -> (Vec<InDel>, Vec<(usize, usize)>) {
  let breakpoints = collect_breakpoints(
    parent_gaps
      .iter()
      .copied()
      .chain(node_gaps.iter().copied())
      .chain(node_non_char.iter().copied())
      .chain(variable_indel.iter().copied()),
  );

  let mut new_node_gaps: Vec<(usize, usize)> = Vec::new();
  let mut deletions: Vec<InDel> = Vec::new();
  let mut insertions: Vec<InDel> = Vec::new();

  for (lo, hi) in breakpoints.iter().copied().tuple_windows() {
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
      if let Some(last) = new_node_gaps.last_mut().filter(|l| l.1 == lo) {
        last.1 = hi;
      } else {
        new_node_gaps.push((lo, hi));
      }
    } else if in_parent && !in_node && !in_non_char && !in_vi {
      // Parent has gap, node has sequence: insertion on this edge.
      if let Some(last) = insertions.last_mut().filter(|i| i.range.1 == lo) {
        last.range.1 = hi;
        last.seq.extend(node_sequence[lo..hi].iter().copied());
      } else {
        insertions.push(InDel::ins((lo, hi), &node_sequence[lo..hi]));
      }
    } else if in_parent {
      // Parent has gap; node is gapped, non_char, or variable: inherit parent gap.
      if let Some(last) = new_node_gaps.last_mut().filter(|l| l.1 == lo) {
        last.1 = hi;
      } else {
        new_node_gaps.push((lo, hi));
      }
    }
  }

  (deletions.into_iter().chain(insertions).collect(), new_node_gaps)
}
