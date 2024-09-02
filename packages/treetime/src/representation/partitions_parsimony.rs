use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::get_common_length;
use crate::io::fasta::FastaRecord;
use eyre::Report;

#[derive(Clone, Debug)]
pub struct PartitionParsimonyWithAln {
  pub alphabet: Alphabet,
  pub aln: Vec<FastaRecord>,
  pub length: usize,
}

impl PartitionParsimonyWithAln {
  pub fn new(alphabet: Alphabet, aln: Vec<FastaRecord>) -> Result<Self, Report> {
    let length = get_common_length(&aln)?;
    Ok(Self { alphabet, aln, length })
  }
}

#[derive(Clone, Debug)]
pub struct PartitionParsimony {
  pub alphabet: Alphabet,
  pub length: usize,
}

impl PartitionParsimony {
  pub fn new(alphabet: Alphabet, length: usize) -> Result<Self, Report> {
    Ok(Self { alphabet, length })
  }
}

impl From<PartitionParsimonyWithAln> for PartitionParsimony {
  fn from(item: PartitionParsimonyWithAln) -> Self {
    Self {
      alphabet: item.alphabet,
      length: item.length,
    }
  }
}
