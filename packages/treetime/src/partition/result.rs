use crate::alphabet::alphabet::Alphabet;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::partition::resolve::PartitionLayout;
use crate::partition::split::AlignmentInfo;
use crate::seq::mutation::{Mutation, SequenceMutations, edge_state_change_counts};
use eyre::Report;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::key_map::EdgeMap;
use treetime_primitives::{LogLh, Seq};
use treetime_utils::make_internal_error;

#[derive(Clone, Debug)]
pub struct PartitionResult {
  pub number: usize,
  pub name: Option<String>,
  pub alignment: usize,
  pub sites: usize,
  pub model_name: GtrModelName,
  pub gtr: Option<GTR>,
  pub log_lh: Option<LogLh>,
}

impl PartitionResult {
  pub(crate) fn new(
    layout: &PartitionLayout,
    model_name: GtrModelName,
    gtr: Option<GTR>,
    log_lh: Option<LogLh>,
  ) -> Self {
    Self {
      number: layout.number,
      name: layout.name.clone(),
      alignment: layout.alignment,
      sites: layout.columns.len(),
      model_name,
      gtr,
      log_lh,
    }
  }
}

pub(crate) fn alignment_results(
  alignments: Vec<AlignmentInfo>,
  outcomes: Vec<AlignmentOutcome>,
) -> Result<Vec<AlignmentResult>, Report> {
  let mut by_alignment = BTreeMap::new();
  for outcome in outcomes {
    let alignment = outcome.alignment;
    if by_alignment.insert(alignment, outcome).is_some() {
      return make_internal_error!("Alignment {alignment} has two results");
    }
  }
  alignments
    .into_iter()
    .map(|info| {
      let Some(AlignmentOutcome {
        indel_rate, sequences, ..
      }) = by_alignment.remove(&info.number)
      else {
        return make_internal_error!("Alignment {} has no result", info.number);
      };
      let AlignmentInfo {
        number,
        path,
        alphabet,
        length,
        mask,
      } = info;
      Ok(AlignmentResult {
        number,
        path,
        alphabet,
        length,
        mask,
        indel_rate,
        sequences,
      })
    })
    .collect()
}

#[derive(Clone, Debug)]
pub struct AlignmentResult {
  pub number: usize,
  pub path: PathBuf,
  pub alphabet: Alphabet,
  pub length: usize,
  pub mask: Vec<bool>,
  pub indel_rate: Option<f64>,
  pub sequences: Option<AlignmentSequences>,
}

pub(crate) struct AlignmentOutcome {
  pub(crate) alignment: usize,
  pub(crate) indel_rate: Option<f64>,
  pub(crate) sequences: Option<AlignmentSequences>,
}

#[derive(Clone, Debug)]
pub struct AlignmentSequences {
  pub root_sequence: Seq,
  pub edge_mutations: EdgeMap<Vec<Mutation>>,
  pub edge_mutation_counts: EdgeMap<usize>,
}

impl AlignmentSequences {
  pub(crate) fn new(mutations: SequenceMutations, alphabet: &Alphabet) -> Result<Self, Report> {
    let SequenceMutations {
      root_sequence,
      edge_mutations,
    } = mutations;
    let edge_mutation_counts = edge_state_change_counts(&edge_mutations, alphabet)?;
    Ok(Self {
      root_sequence,
      edge_mutations,
      edge_mutation_counts,
    })
  }
}
