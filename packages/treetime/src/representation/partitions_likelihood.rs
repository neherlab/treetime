use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::io::fasta::FastaRecord;
use crate::port::fitch::get_common_length;
use crate::representation::graph_sparse::VarPos;
use crate::representation::partitions_parsimony::PartitionParsimony;
use crate::seq::range::RangeCollection;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

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

impl From<PartitionLikelihoodWithAln> for PartitionLikelihood {
  fn from(item: PartitionLikelihoodWithAln) -> Self {
    Self {
      gtr: item.gtr,
      alphabet: item.alphabet,
      length: item.length,
    }
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

  pub fn from_parsimony(gtr: GTR, partition: PartitionParsimony) -> Self {
    Self::new(gtr, partition.alphabet, partition.length)
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SeqInfoLikelihood {
  pub unknown: RangeCollection,
  pub gaps: RangeCollection,
  pub non_char: RangeCollection, // any position that does not evolve according to the substitution model, i.e. gap or N
  pub variable: BTreeMap<usize, VarPos>, // probability vector for each variable position collecting information from children
  pub log_lh: f64,
  pub fixed: BTreeMap<String, Array1<f64>>, // probability vector for the state of fixed positions based on information from children
  pub fixed_composition: BTreeMap<String, usize>, // this is important for likelihood calculations
}
