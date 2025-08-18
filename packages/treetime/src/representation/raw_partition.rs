use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::io::fasta::FastaRecord;
use crate::representation::partitions_likelihood::PartitionLikelihoodWithAln;
use crate::representation::partitions_parsimony::PartitionParsimonyWithAln;
use enum_extract_macro::EnumExtract;
use eyre::Report;

#[derive(Clone, Debug, EnumExtract)]
pub enum RawPartition {
  Parsimony(PartitionParsimonyWithAln),
  Likelihood(PartitionLikelihoodWithAln),
}

impl RawPartition {
  pub fn parsimony(alphabet: Alphabet, aln: Vec<FastaRecord>) -> Result<Self, Report> {
    Ok(RawPartition::Parsimony(PartitionParsimonyWithAln::new(alphabet, aln)?))
  }

  pub fn likelihood(gtr: GTR, alphabet: Alphabet, aln: Vec<FastaRecord>) -> Result<Self, Report> {
    Ok(RawPartition::Likelihood(PartitionLikelihoodWithAln::new(
      gtr, alphabet, aln,
    )?))
  }
}
