use crate::ancestral::params::AncestralParams;
use crate::ancestral::plan::{ReconstructionOptions, reconstruct_partitions, resolve_plan};
use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::partition::result::{PartitionResult, PartitionSequences, split_sources};
use crate::progress::{LogSink, StageSink};
use crate::seq::alignment::AncestralInput;
use crate::seq::sink::{SeqSink, SeqTrack};
use treetime_graph::graph::Graph;
use treetime_utils::make_report;
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
  if input.sources.is_empty() {
    return Err(OperationError::InvalidInput(make_report!(
      "Ancestral reconstruction requires at least one alignment"
    )));
  }
  let options = ReconstructionOptions::new(params.impute_missing_data, params.sample_from_profile);
  let mut rng = get_random_number_generator(params.seed);
  let branch_lengths = branch_lengths_or_zero(&input.branch_lengths());
  let AncestralInput { graph, sources, .. } = input;
  let partitions = reconstruct_partitions(
    &graph,
    &plan,
    split_sources(sources),
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
  let results = partitions
    .into_iter()
    .map(|(source, partition)| {
      let mutations = partition.stream_sequences(
        &graph,
        SeqTrack::Nuc(source.number),
        params.include_leaves,
        params.report_ambiguous,
        seq_sink.as_deref_mut(),
      )?;
      let alphabet = partition.alphabet().clone();
      let sequences = PartitionSequences::new(mutations, &alphabet).map_err(OperationError::InferenceFailed)?;
      Ok(PartitionResult {
        number: source.number,
        path: source.path,
        length: source.length,
        mask: source.mask,
        model_name: params.model,
        gtr: partition.gtr().cloned(),
        indel_rate: None,
        log_lh: partition.log_lh(&graph).map_err(OperationError::classify)?,
        sequences: Some(sequences),
        alphabet,
      })
    })
    .collect::<Result<Vec<_>, OperationError>>()?;
  Ok(AncestralOutput { results, graph })
}

#[derive(Debug)]
pub struct AncestralOutput {
  pub results: Vec<PartitionResult>,
  pub graph: Graph,
}
