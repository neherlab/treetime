use crate::alphabet::alphabet::Alphabet;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::partition::resolve::PartitionLayout;
use crate::partition::split::AlignmentInfo;
use crate::seq::mutation::{Mutation, SequenceMutations, edge_state_change_counts};
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::edge::GraphEdgeKey;
use treetime_primitives::{LogLh, Seq};
use treetime_utils::make_internal_error;

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
  outcomes: Vec<PartitionOutcome>,
) -> Result<Vec<AlignmentResult>, Report> {
  let mut by_alignment = outcomes.into_iter().into_group_map_by(|outcome| outcome.alignment);
  alignments
    .into_iter()
    .map(|info| {
      let outcomes = by_alignment.remove(&info.number).unwrap_or_default();
      let indel_rate = outcomes
        .iter()
        .filter_map(|outcome| outcome.indel_rate)
        .reduce(|total, rate| total + rate);
      let sequences = match <[PartitionOutcome; 1]>::try_from(outcomes) {
        Ok([outcome]) => outcome.sequences,
        Err(outcomes) if outcomes.iter().all(|outcome| outcome.sequences.is_none()) => None,
        Err(outcomes) => {
          return make_internal_error!(
            "Alignment {} has the sequences of {} partitions, but sequences of a split alignment cannot be joined",
            info.number,
            outcomes.len()
          );
        },
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

pub(crate) struct PartitionOutcome {
  pub(crate) alignment: usize,
  pub(crate) indel_rate: Option<f64>,
  pub(crate) sequences: Option<AlignmentSequences>,
}

#[derive(Clone, Debug)]
pub struct AlignmentSequences {
  pub root_sequence: Seq,
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  pub edge_mutation_counts: BTreeMap<GraphEdgeKey, usize>,
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
