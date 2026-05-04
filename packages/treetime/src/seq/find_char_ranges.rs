use treetime_primitives::{AsciiChar, Seq};

// Finds maximal contiguous ranges in the sequence where the predicate returns true for every character.
#[inline]
pub fn find_letter_ranges_by<F>(seq: &Seq, pred: F) -> Vec<(usize, usize)>
where
  F: Fn(AsciiChar) -> bool + Copy,
{
  let mut result = Vec::with_capacity(31);
  let mut i = 0;
  while i < seq.len() {
    if pred(seq[i]) {
      let start = i;
      i += 1;
      while i < seq.len() && pred(seq[i]) {
        i += 1;
      }
      result.push((start, i));
    } else {
      i += 1;
    }
  }
  result
}

/// Finds contiguous ranges (segments) consisting of a given letter in the sequence.
#[inline]
pub fn find_letter_ranges(seq: &Seq, letter: AsciiChar) -> Vec<(usize, usize)> {
  find_letter_ranges_by(seq, |candidate: AsciiChar| candidate == letter)
}
