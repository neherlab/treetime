// Finds contiguous ranges (segments) in the sequence, such that for every character inside every range,
// the predicate function returns true and every range contains only the same letter.
//
// The predicate is a function that takes a character and returns boolean.
//
// For example if predicate returns `true` for characters A and C, this function will find ranges `AAAA` and `CCCCC`,
// but not `ZZZ` or `ACCCAC`.
pub fn find_letter_ranges_by(seq: &[char], pred: impl Fn(char) -> bool) -> Vec<(usize, usize)> {
  let len = seq.len();

  let mut result = vec![];
  let mut i = 0_usize;
  let mut start = 0_usize;
  let mut found_maybe = Option::<char>::default();
  while i < len {
    let letter = seq[i];

    // Find beginning of a range
    if pred(letter) {
      start = i;
      found_maybe = Some(letter);
    }

    match found_maybe {
      // If there's a current range we are working on (for which we found a `start`), extend it
      Some(found) => {
        // Rewind forward until we find a mismatch
        while i < len && seq[i] == found {
          i += 1;
        }

        // We found the end of the current range, so now it's complete
        let end = i;

        // Remember the range
        result.push((start, end));

        found_maybe = None;
      }
      None => {
        if i < len {
          i += 1;
        }
      }
    }
  }
  result
}

/// Finds contiguous ranges (segments) consisting of a given letter in the sequence.
pub fn find_letter_ranges(seq: &[char], letter: char) -> Vec<(usize, usize)> {
  find_letter_ranges_by(seq, |candidate| candidate == letter)
}

pub fn find_ambiguous_ranges(seq: &[char]) -> Vec<(usize, usize)> {
  find_letter_ranges(seq, 'N')
}

pub fn find_gap_ranges(seq: &[char]) -> Vec<(usize, usize)> {
  find_letter_ranges(seq, '-')
}
