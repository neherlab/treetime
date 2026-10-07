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
use treetime_graph::graph::Graph;
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::sync::random::get_random_number_generator;

pub fn run(
  params: &AncestralParams,
  input: AncestralInput,
  mut seq_sink: Option<&mut dyn SeqSink>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<AncestralOutput, OperationError> {
  let plan = resolve_plan(params)?;
  let options = ReconstructionOptions::new(params.impute_missing_data, params.sample_from_profile);
  let mut rng = get_random_number_generator(params.seed);
  let branch_lengths = branch_lengths_or_zero(&input.branch_lengths());
  let AncestralInput {
    graph,
    nodes,
    alphabet,
    mask,
    ..
  } = input;
  let ambiguous_char = alphabet.unknown();
  let partition = reconstruct_partition(
    &graph,
    &plan,
    alphabet,
    nodes,
    &branch_lengths,
    &options,
    &mut rng,
    cancel,
    stages,
    log,
  )
  .map_err(OperationError::classify)?;
  if let Some(sink) = seq_sink.as_deref_mut() {
    sink.on_topology(&graph).map_err(OperationError::SinkFailed)?;
  }
  let SequenceMutations {
    root_sequence,
    edge_mutations,
  } = partition.stream_sequences(
    &graph,
    &MutationTrack::Nucleotide,
    params.include_leaves,
    params.report_ambiguous,
    seq_sink,
  )?;
  Ok(AncestralOutput {
    gtr: partition.gtr().cloned(),
    model_name: params.model,
    method: params.method,
    mask,
    sequence_length: root_sequence.len(),
    ambiguous_char,
    root_sequence,
    edge_mutations,
    graph,
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
  pub graph: Graph,
}
