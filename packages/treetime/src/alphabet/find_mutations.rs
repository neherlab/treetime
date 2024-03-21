use crate::alphabet::sequence_data::SequenceData;
use crate::commands::ancestral::anc_graph::AncestralGraph;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::graph::GraphNodeForward;
use ndarray::{Array1, CowArray, Zip};
use num_traits::{One, Zero};
use std::cmp::Ordering;
use std::fmt::{Display, Formatter};

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct Mutation {
  pub reff: char,
  pub pos: usize,
  pub qry: char,
}

impl Display for Mutation {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}{}{}", self.reff, self.pos, self.qry)
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

  pub const fn is_valid(&self) -> bool {
    let Mutation { reff, pos, qry } = self;
    *reff != '-' && *reff != 'N' && *qry != '-' && *qry != 'N'
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

/// Finds mutations on a graph and attach them to node payloads.
pub fn find_graph_mutations(graph: &mut AncestralGraph, sequence_data: &SequenceData) {
  graph.par_iter_breadth_first_forward(
    |GraphNodeForward {
       is_root,
       is_leaf,
       key,
       payload: mut node,
       parents,
     }| {
      if is_root {
        return GraphTraversalContinuation::Continue;
      }

      if parents.len() > 1 {
        unimplemented!("Multiple parent nodes are not supported yet");
      }

      let (parent, _) = &parents[0];

      let parent_seq = sequence_data.get_full(&parent.read().name).unwrap();
      let this_seq = sequence_data.get_full(&node.name).unwrap();
      node.mutations = find_mutations(&parent_seq, &this_seq, &node.mask);

      GraphTraversalContinuation::Continue
    },
  );
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
