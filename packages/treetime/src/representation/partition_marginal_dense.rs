use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;

#[derive(Clone, Debug)]
pub struct PartitionMarginalDense {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
}
