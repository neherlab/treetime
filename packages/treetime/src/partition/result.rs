use crate::alphabet::alphabet::Alphabet;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::seq::alignment::{NodeSeqInput, SequenceSource};
use crate::seq::mutation::{Mutation, SequenceMutations, edge_state_change_counts};
use eyre::Report;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{LogLh, Seq};

#[derive(Clone, Debug)]
pub struct PartitionResult {
  pub number: usize,
  pub path: PathBuf,
  pub alphabet: Alphabet,
  pub length: usize,
  pub mask: Vec<bool>,
  pub model_name: GtrModelName,
  pub gtr: Option<GTR>,
  pub indel_rate: Option<f64>,
  pub log_lh: Option<LogLh>,
  pub sequences: Option<PartitionSequences>,
}

#[derive(Clone, Debug)]
pub struct PartitionSequences {
  pub root_sequence: Seq,
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  pub edge_mutation_counts: BTreeMap<GraphEdgeKey, usize>,
}

impl PartitionSequences {
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

pub(crate) fn split_sources(sources: Vec<SequenceSource>) -> Vec<SplitSource> {
  sources
    .into_iter()
    .enumerate()
    .map(|(index, source)| split_source(index + 1, source))
    .collect()
}

pub(crate) fn split_source(number: usize, source: SequenceSource) -> SplitSource {
  let SequenceSource {
    path,
    alphabet,
    length,
    nodes,
    mask,
  } = source;
  SplitSource {
    source: PartitionSource {
      number,
      path,
      length,
      mask,
    },
    alphabet,
    nodes,
  }
}

pub(crate) struct SplitSource {
  pub(crate) source: PartitionSource,
  pub(crate) alphabet: Alphabet,
  pub(crate) nodes: BTreeMap<GraphNodeKey, NodeSeqInput>,
}

pub(crate) struct PartitionSource {
  pub(crate) number: usize,
  pub(crate) path: PathBuf,
  pub(crate) length: usize,
  pub(crate) mask: Vec<bool>,
}
