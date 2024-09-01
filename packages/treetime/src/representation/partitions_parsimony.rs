use crate::alphabet::alphabet::Alphabet;
use crate::io::fasta::FastaRecord;
use crate::port::fitch::get_common_length;
use crate::representation::graph_sparse::VarPos;
use crate::seq::range::RangeCollection;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

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

impl From<PartitionParsimonyWithAln> for PartitionParsimony {
  fn from(item: PartitionParsimonyWithAln) -> Self {
    Self {
      alphabet: item.alphabet,
      length: item.length,
    }
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

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SeqInfoParsimony {
  pub unknown: RangeCollection,
  pub gaps: RangeCollection,
  pub non_char: RangeCollection, // any position that does not evolve according to the substitution model, i.e. gap or N
  pub variable: BTreeMap<usize, VarPos>, // 0/1 vector for each variable position collecting information from children
  pub fixed_composition: BTreeMap<String, usize>, // this is important for likelihood calculations
}
