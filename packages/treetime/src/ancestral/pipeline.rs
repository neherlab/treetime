use crate::ancestral::params::{AncestralParams, MethodAncestral};
use crate::ancestral::plan::{ReconstructionOptions, reconstruct_partition, resolve_plan};
use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::progress::{LogSink, StageSink};
use crate::seq::alignment::AncestralInput;
use crate::seq::mutation::{Mutation, MutationTrack, SequenceMutations};
use crate::seq::sink::SeqSink;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::sync::random::get_random_number_generator;

pub fn run(
  params: &AncestralParams,
  input: &AncestralInput,
  seq_sink: &mut dyn SeqSink,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<AncestralOutput, OperationError> {
  let plan = resolve_plan(params)?;
  let options = ReconstructionOptions::new(params.impute_missing_data, params.sample_from_profile);
  let mut rng = get_random_number_generator(params.seed);
  let branch_lengths = branch_lengths_or_zero(&input.branch_lengths());
  let partition = reconstruct_partition(
    &input.graph,
    &plan,
    0,
    input.alphabet.clone(),
    &input.nodes,
    &branch_lengths,
    &options,
    &mut rng,
    cancel,
    stages,
    log,
  )
  .map_err(OperationError::classify)?;
  seq_sink.on_topology(&input.graph).map_err(OperationError::SinkFailed)?;
  let SequenceMutations {
    root_sequence,
    edge_mutations,
  } = partition.stream_sequences(
    &input.graph,
    &MutationTrack::Nucleotide,
    params.include_leaves,
    Some(seq_sink),
  )?;
  Ok(AncestralOutput {
    gtr: partition.gtr().cloned(),
    model_name: params.model,
    method: params.method,
    mask: input.mask.clone(),
    sequence_length: root_sequence.len(),
    ambiguous_char: input.alphabet.unknown(),
    root_sequence,
    edge_mutations,
  })
}

#[derive(Debug)]
pub struct AncestralOutput {
  pub gtr: Option<GTR>,
  pub model_name: GtrModelName,
  pub method: MethodAncestral,
  pub mask: Vec<bool>,
  pub sequence_length: usize,
  pub ambiguous_char: AsciiChar,
  pub root_sequence: Seq,
  pub edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}
