use crate::alphabet::alphabet::Alphabet;
use treetime_primitives::AlignmentRecord;

pub fn create_mask(aln: &[AlignmentRecord], alignment_length: usize, alphabet: &Alphabet) -> Vec<bool> {
  let ambiguous = alphabet.unknown();
  let gap = alphabet.gap();

  (0..alignment_length)
    .map(|pos| {
      !aln.iter().any(|record| {
        record
          .seq
          .as_slice()
          .get(pos)
          .is_some_and(|&state| state != ambiguous && state != gap)
      })
    })
    .collect()
}

pub fn mask_to_string(mask: &[bool]) -> String {
  mask.iter().map(|&m| if m { '1' } else { '0' }).collect()
}
