use crate::make_internal_error;
use eyre::Report;
use maplit::btreeset;
use std::collections::{BTreeMap, BTreeSet};

pub fn find_mixed_sites(seq: &[char]) -> (Vec<MixedSite>, BTreeMap<usize, BTreeSet<char>>) {
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

  pub fn disambiguate(&self) -> Result<BTreeSet<char>, Report> {
    match &self.nuc {
      'R' => Ok(btreeset! {'A', 'G'}),
      'Y' => Ok(btreeset! {'C', 'T'}),
      'S' => Ok(btreeset! {'G', 'C'}),
      'W' => Ok(btreeset! {'A', 'T'}),
      'K' => Ok(btreeset! {'G', 'T'}),
      'M' => Ok(btreeset! {'A', 'C'}),
      'B' => Ok(btreeset! {'C', 'G', 'T'}),
      'D' => Ok(btreeset! {'A', 'G', 'T'}),
      'H' => Ok(btreeset! {'A', 'C', 'T'}),
      'V' => Ok(btreeset! {'A', 'C', 'G'}),
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
      (3, btreeset! {'A', 'G', 'T'}),
      (9, btreeset! {'C', 'T'}),
      (10, btreeset! {'A', 'G'}),
      (14, btreeset! {'A', 'C', 'T'}),
    ]);

    let expected = (mixed, disambiguated);

    assert_eq!(expected, actual);
    Ok(())
  }
}
