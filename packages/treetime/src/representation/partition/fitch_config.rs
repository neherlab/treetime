use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::get_common_length;
use eyre::Report;
use treetime_io::fasta::FastaRecord;

#[derive(Clone, Debug)]
pub struct PartitionFitchConfigWithAln {
  pub alphabet: Alphabet,
  pub aln: Vec<FastaRecord>,
  pub length: usize,
}

impl PartitionFitchConfigWithAln {
  pub fn new(alphabet: Alphabet, aln: Vec<FastaRecord>) -> Result<Self, Report> {
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
