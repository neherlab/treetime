use log::debug;
use serde::{Deserialize, Serialize};
use std::fmt;
use treetime_primitives::Seq;

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
    debug!(
      "InDel::new: kind={}, range={:?}, seq={}",
      if deletion { "del" } else { "ins" },
      range,
      seq.as_str(),
    );
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
