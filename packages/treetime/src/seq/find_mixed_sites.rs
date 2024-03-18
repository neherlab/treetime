use crate::make_internal_error;
use eyre::Report;
use std::collections::BTreeMap;

pub fn find_mixed_sites(seq: &[char]) -> (Vec<MixedSite>, BTreeMap<usize, Vec<char>>) {
  let mut mixed_positions = Vec::new();

  for (pos, nuc) in seq.iter().enumerate() {
    if !"ACGTN-".contains(*nuc) {
      mixed_positions.push(MixedSite::new(pos, *nuc));
    }
  }

  let non_consensus = mixed_positions
    .iter()
    .map(|ms| (ms.pos, ms.disambiguate().unwrap()))
    .collect();

  (mixed_positions, non_consensus)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MixedSite {
  pub pos: usize,
  pub nuc: char,
}

impl MixedSite {
  pub fn new(pos: usize, nuc: char) -> Self {
    Self { pos, nuc }
  }

  pub fn disambiguate(&self) -> Result<Vec<char>, Report> {
    match &self.nuc {
      'R' => Ok(vec!['A', 'G']),
      'Y' => Ok(vec!['C', 'T']),
      'S' => Ok(vec!['G', 'C']),
      'W' => Ok(vec!['A', 'T']),
      'K' => Ok(vec!['G', 'T']),
      'M' => Ok(vec!['A', 'C']),
      'B' => Ok(vec!['C', 'G', 'T']),
      'D' => Ok(vec!['A', 'G', 'T']),
      'H' => Ok(vec!['A', 'C', 'T']),
      'V' => Ok(vec!['A', 'C', 'G']),
      _ => make_internal_error!("Unknown ambiguous nucleotide: '{}'", &self.nuc),
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn test_find_mixed_sites() -> Result<(), Report> {
    //                              012345678901234
    let actual = find_mixed_sites(&"TCGDCNGTGYRTTGH".chars().collect_vec());
    //                                 D     YR   H

    let mixed = vec![
      MixedSite::new(3, 'D'),
      MixedSite::new(9, 'Y'),
      MixedSite::new(10, 'R'),
      MixedSite::new(14, 'H'),
    ];

    let disambiguated = BTreeMap::from([
      (3, vec!['A', 'G', 'T']),
      (9, vec!['C', 'T']),
      (10, vec!['A', 'G']),
      (14, vec!['A', 'C', 'T']),
    ]);

    let expected = (mixed, disambiguated);

    assert_eq!(expected, actual);
    Ok(())
  }
}
