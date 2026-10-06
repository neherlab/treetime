use crate::alphabet::alphabet::Alphabet;
use itertools::Itertools;
use treetime_primitives::Seq;

pub fn create_mask<'a>(
  seqs: impl IntoIterator<Item = &'a Seq>,
  alignment_length: usize,
  alphabet: &Alphabet,
) -> Vec<bool> {
  let ambiguous = alphabet.unknown();
  let gap = alphabet.gap();
  let seqs = seqs.into_iter().collect_vec();

  (0..alignment_length)
    .map(|pos| {
      !seqs.iter().any(|seq| {
        seq
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
