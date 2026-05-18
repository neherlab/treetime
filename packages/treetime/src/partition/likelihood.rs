use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::partition::fitch_config::PartitionFitchConfig;
use crate::seq::alignment::get_common_length;
use eyre::Report;
use treetime_io::fasta::FastaRecord;

#[derive(Clone, Debug)]
pub struct PartitionLikelihoodWithAln {
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub aln: Vec<FastaRecord>,
  pub length: usize,
}

impl PartitionLikelihoodWithAln {
  pub fn new(gtr: GTR, alphabet: Alphabet, aln: Vec<FastaRecord>) -> Result<Self, Report> {
    let length = get_common_length(&aln)?;
    Ok(Self {
      gtr,
      alphabet,
      aln,
      length,
    })
  }
}

#[derive(Clone, Debug)]
pub struct PartitionLikelihood {
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
}

impl PartitionLikelihood {
  pub fn new(gtr: GTR, alphabet: Alphabet, length: usize) -> Self {
    Self { gtr, alphabet, length }
  }

  pub fn from_fitch(gtr: GTR, partition: PartitionFitchConfig) -> Self {
    Self::new(gtr, partition.alphabet, partition.length)
  }
}

impl From<PartitionLikelihoodWithAln> for PartitionLikelihood {
  fn from(item: PartitionLikelihoodWithAln) -> Self {
    Self {
      gtr: item.gtr,
      alphabet: item.alphabet,
      length: item.length,
    }
  }
}
