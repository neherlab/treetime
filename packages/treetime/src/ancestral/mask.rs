use crate::alphabet::alphabet::Alphabet;
use treetime_io::fasta::FastaRecord;

pub fn create_mask(aln: &[FastaRecord], alignment_length: usize, alphabet: &Alphabet) -> Vec<bool> {
  let ambiguous = alphabet.unknown();
  let gap = alphabet.gap();

  let mut mask = vec![true; alignment_length];

  for record in aln {
    for (pos, &state) in record.seq.iter().enumerate() {
      if pos < alignment_length && state != ambiguous && state != gap {
        mask[pos] = false;
      }
    }
  }

  mask
}

pub fn mask_to_string(mask: &[bool]) -> String {
  mask.iter().map(|&m| if m { '1' } else { '0' }).collect()
}
