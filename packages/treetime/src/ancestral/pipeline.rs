use crate::ancestral::params::AncestralParams;
use crate::ancestral::plan::{ReconstructedPartition, ReconstructionOptions, reconstruct_partitions, resolve_plan};
use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::partition::result::{
  AlignmentResult, AlignmentSequences, PartitionOutcome, PartitionResult, alignment_results,
};
use crate::progress::{LogSink, StageSink};
use crate::seq::alignment::AncestralInput;
use crate::seq::sink::{SeqStream, SeqTrack};
use treetime_graph::graph::Graph;
use treetime_utils::make_report;
use treetime_utils::sync::random::get_random_number_generator;

pub fn run(
  params: &AncestralParams,
  input: AncestralInput,
  mut seq_stream: Option<SeqStream<'_>>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<AncestralOutput, OperationError> {
  let method = resolve_plan(params)?;
  if input.sequences.is_empty() {
    return Err(OperationError::InvalidInput(make_report!(
      "Ancestral reconstruction requires at least one alignment"
    )));
  }
  let options = ReconstructionOptions::new(params.impute_missing_data, params.sample_from_profile);
  let mut rng = get_random_number_generator(params.seed);
  let branch_lengths = branch_lengths_or_zero(&input.branch_lengths());
  let AncestralInput { graph, sequences, .. } = input;
  let split = sequences.split().map_err(OperationError::classify)?;
  let alignments = split.alignments.clone();
  let reconstructed = reconstruct_partitions(
    &graph,
    method,
    split,
    &branch_lengths,
    &options,
    &mut rng,
    cancel,
    stages,
    log,
  )
  .map_err(OperationError::classify)?;
  let (partitions, outcomes): (Vec<_>, Vec<_>) = reconstructed
    .into_iter()
    .map(
      |ReconstructedPartition {
         layout,
         model_name,
         partition,
       }| {
        let mutations = partition.stream_sequences(
          &graph,
          SeqTrack::Nuc(layout.alignment),
          params.include_leaves,
          params.report_ambiguous,
          seq_stream.as_mut().map(SeqStream::reborrow),
        )?;
        let sequences =
          AlignmentSequences::new(mutations, partition.alphabet()).map_err(OperationError::InferenceFailed)?;
        let log_lh = partition.log_lh(&graph).map_err(OperationError::classify)?;
        let result = PartitionResult::new(&layout, model_name, partition.gtr().cloned(), log_lh);
        let outcome = PartitionOutcome {
          alignment: layout.alignment,
          indel_rate: None,
          sequences: Some(sequences),
        };
        Ok((result, outcome))
      },
    )
    .collect::<Result<Vec<_>, OperationError>>()?
    .into_iter()
    .unzip();
  let alignments = alignment_results(alignments, outcomes).map_err(OperationError::InferenceFailed)?;
  Ok(AncestralOutput {
    alignments,
    partitions,
    graph,
  })
}

#[derive(Debug)]
pub struct AncestralOutput {
  pub alignments: Vec<AlignmentResult>,
  pub partitions: Vec<PartitionResult>,
  pub graph: Graph,
}
