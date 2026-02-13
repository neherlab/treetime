use treetime_primitives::{AsciiChar, Seq};

// Finds contiguous ranges (segments) in the sequence, such that for every character inside every range,
// the predicate function returns true and every range contains only the same letter.
//
// The predicate is a function that takes a character and returns boolean.
//
// For example if predicate returns `true` for characters A and C, this function will find ranges `AAAA`, `CCCCC`, `ACCCACAAA`
// but not `ZZZZZ`.
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

pub mod old {
  #![allow(dead_code)]

  /// OLD, SLOWER VERSION
  pub fn find_letter_ranges_by(seq: &[char], pred: impl Fn(char) -> bool) -> Vec<(usize, usize)> {
    let len = seq.len();

    let mut result = vec![];
    let mut i = 0_usize;
    let mut start = 0_usize;
    let mut is_inside_range = false;
    while i < len {
      let letter = seq[i];

      // Find beginning of a range
      if pred(letter) {
        start = i;
        is_inside_range = true;
      }

      if is_inside_range {
        // Rewind forward until we find a mismatch
        while i < len && pred(seq[i]) {
          i += 1;
        }

        // We found the end of the current range, so now it's complete
        let end = i;

        // Remember the range
        result.push((start, end));

        is_inside_range = false;
      } else if i < len {
        i += 1;
      }
    }
    result
  }
}

/// Finds contiguous ranges (segments) consisting of a given letter in the sequence.
#[inline]
pub fn find_letter_ranges(seq: &Seq, letter: impl Into<AsciiChar> + Copy) -> Vec<(usize, usize)> {
  find_letter_ranges_by(seq, |candidate: AsciiChar| candidate == letter.into())
}
