use crate::alphabet::alphabet::Alphabet;
use crate::seq::alignment::get_common_length;
use eyre::Report;
use treetime_primitives::AlignmentRecord;

#[derive(Clone, Debug)]
pub struct PartitionFitchConfigWithAln {
  pub alphabet: Alphabet,
  pub aln: Vec<AlignmentRecord>,
  pub length: usize,
}

impl PartitionFitchConfigWithAln {
  pub fn new(alphabet: Alphabet, aln: Vec<AlignmentRecord>) -> Result<Self, Report> {
    let length = get_common_length(&aln)?;
    Ok(Self { alphabet, aln, length })
  }
}

#[derive(Clone, Debug)]
pub struct PartitionFitchConfig {
  pub alphabet: Alphabet,
  pub length: usize,
}

impl PartitionFitchConfig {
  pub fn new(alphabet: Alphabet, length: usize) -> Self {
    Self { alphabet, length }
  }
}

impl From<PartitionFitchConfigWithAln> for PartitionFitchConfig {
  fn from(item: PartitionFitchConfigWithAln) -> Self {
    Self {
      alphabet: item.alphabet,
      length: item.length,
    }
  }
}
