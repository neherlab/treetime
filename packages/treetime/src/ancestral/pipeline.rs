use crate::ancestral::params::AncestralParams;
use crate::ancestral::plan::{ReconstructedPartition, ReconstructionOptions, reconstruct_partitions, resolve_plan};
use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::partition::alignment_sequences::{StreamOptions, stream_alignment_sequences};
use crate::partition::assembly::alignment_groups;
use crate::partition::result::{AlignmentOutcome, AlignmentResult, PartitionResult, alignment_results};
use crate::progress::{LogSink, StageSink};
use crate::seq::alignment::AncestralInput;
use crate::seq::sink::SeqStream;
use itertools::Itertools;
use treetime_graph::graph::Graph;
use treetime_utils::make_report;
use treetime_utils::sync::random::get_random_number_generator;

pub fn run(
  params: &AncestralParams,
  input: AncestralInput,
  seq_stream: Option<SeqStream<'_>>,
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
  let (outcomes, partitions) = ancestral_results(params, &graph, reconstructed, seq_stream)?;
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

fn ancestral_results(
  params: &AncestralParams,
  graph: &Graph,
  reconstructed: Vec<ReconstructedPartition>,
  mut seq_stream: Option<SeqStream<'_>>,
) -> Result<(Vec<AlignmentOutcome>, Vec<PartitionResult>), OperationError> {
  let groups = alignment_groups(
    reconstructed
      .iter()
      .map(|partition| (Some(&partition.layout), partition.layout.columns.len())),
  )
  .map_err(OperationError::InferenceFailed)?;
  let stream_options = StreamOptions {
    include_leaves: params.include_leaves,
    report_unknown: params.report_ambiguous,
  };
  let reconstructions = reconstructed
    .iter()
    .map(|reconstructed| &reconstructed.partition)
    .collect_vec();
  let outcomes = groups
    .iter()
    .map(|group| {
      let sequences = stream_alignment_sequences(
        graph,
        group,
        &group.members(&reconstructions),
        stream_options,
        seq_stream.as_mut().map(SeqStream::reborrow),
      )?;
      Ok(AlignmentOutcome {
        alignment: group.alignment(),
        indel_rate: None,
        sequences: Some(sequences),
      })
    })
    .collect::<Result<Vec<_>, OperationError>>()?;
  let partitions = reconstructed
    .into_iter()
    .map(
      |ReconstructedPartition {
         layout,
         model_name,
         partition,
       }| {
        let log_lh = partition.log_lh(graph).map_err(OperationError::classify)?;
        Ok(PartitionResult::new(
          &layout,
          model_name,
          partition.gtr().cloned(),
          log_lh,
        ))
      },
    )
    .collect::<Result<Vec<_>, OperationError>>()?;
  Ok((outcomes, partitions))
}
