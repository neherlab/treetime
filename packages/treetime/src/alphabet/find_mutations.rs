use ndarray::{Array1, CowArray, Zip};
use num_traits::{One, Zero};
use std::cmp::Ordering;
use std::fmt::Display;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Mutation {
  reff: char,
  pos: usize,
  qry: char,
}

impl ToString for Mutation {
  fn to_string(&self) -> String {
    format!("{}{}{}", self.reff, self.pos, self.qry)
  }
}

impl Ord for Mutation {
  fn cmp(&self, other: &Self) -> Ordering {
    (self.pos, self.reff, self.qry).cmp(&(other.pos, other.reff, other.qry))
  }
}

impl PartialOrd for Mutation {
  fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
    Some(self.cmp(other))
  }
}

impl Mutation {
  pub const fn new(reff: char, pos: usize, qry: char) -> Self {
    Self { reff, pos, qry }
  }
}

pub fn find_mutations<M: Clone + PartialOrd + Zero + One>(
  ref_seq: &Array1<char>,
  qry_seq: &Array1<char>,
  mask: &Option<Array1<M>>,
) -> Vec<Mutation> {
  assert_eq!(
    ref_seq.len(),
    qry_seq.len(),
    "find_mutations: sequences should have the same length"
  );

  let mask = match mask {
    None => CowArray::from(Array1::<M>::ones(ref_seq.len())),
    Some(mask) => CowArray::from(mask),
  };

  assert_eq!(
    ref_seq.len(),
    mask.len(),
    "find_mutations: mask should have the same length as  reference"
  );

  let indices = Array1::<usize>::from_iter(0..ref_seq.len());
  Zip::from(ref_seq)
    .and(&indices)
    .and(qry_seq)
    .and(mask.view())
    .fold(vec![], |mut acc, reff, pos, qry, mask| {
      if (mask > &M::zero()) && (reff != qry) {
        acc.push(Mutation::new(*reff, *pos, *qry));
      }
      acc
    })
}

#[cfg(test)]
mod tests {
  use super::{find_mutations, Mutation};
  use eyre::Report;
  use ndarray::{array, Array1};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  fn finds_mutations() -> Result<(), Report> {
    //                 0    1    2    3    4    5    6    7
    let reff = array!['C', 'G', 'G', 'T', 'A', 'T', 'C', 'A'];
    let qryy = array!['T', 'G', 'A', 'T', 'A', 'T', 'C', 'G'];
    let expected = vec![
      Mutation::new('C', 0, 'T'),
      Mutation::new('G', 2, 'A'),
      Mutation::new('A', 7, 'G'),
    ];
    let actual = find_mutations(&reff, &qryy, &None::<Array1<i8>>);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rstest]
  fn finds_mutations_with_mask() -> Result<(), Report> {
    //                 0    1    2    3    4    5    6    7
    let reff = array!['C', 'G', 'G', 'T', 'A', 'T', 'C', 'A'];
    let qryy = array!['T', 'G', 'A', 'T', 'A', 'T', 'C', 'G'];
    #[rustfmt::skip]
    let mask = array![ 0 ,  1,   0,   1,   2,   1,   2,   3 ];
    let expected = vec![Mutation::new('A', 7, 'G')];
    let actual = find_mutations(&reff, &qryy, &Some(mask));
    assert_eq!(expected, actual);
    Ok(())
  }
}
