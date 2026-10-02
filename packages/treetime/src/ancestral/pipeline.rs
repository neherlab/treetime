use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::params::AncestralParams;
use crate::ancestral::partition::AncestralPartition;
use crate::ancestral::plan::{ReconstructedPartition, ReconstructionOptions, reconstruct_partition, resolve_plan};
use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::progress::{LogSink, StageSink};
use crate::seq::alignment::AncestralInput;
use serde::Serialize;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::sync::random::get_random_number_generator;

pub fn run(
  params: &AncestralParams,
  input: &AncestralInput,
  alphabet: Alphabet,
  mask: Vec<bool>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<AncestralOutputFull, OperationError> {
  let plan = resolve_plan(params)?;
  let options = ReconstructionOptions::new(
    params.include_leaves,
    params.impute_missing_data,
    params.sample_from_profile,
  );
  let mut rng = get_random_number_generator(params.seed);
  let branch_lengths = branch_lengths_or_zero(&input.branch_lengths());
  let ReconstructedPartition {
    partition,
    emitted_nodes,
  } = reconstruct_partition(
    &input.graph,
    &plan,
    0,
    alphabet,
    &input.nodes,
    &branch_lengths,
    &options,
    &mut rng,
    cancel,
    stages,
    log,
  )?;
  Ok(AncestralOutputFull {
    output: AncestralOutput {
      gtr: partition.gtr().cloned(),
      model_name: params.model,
      mask,
      emitted_nodes,
    },
    partition: Some(partition),
  })
}

pub struct AncestralOutputFull {
  pub output: AncestralOutput,
  pub partition: Option<AncestralPartition>,
}

#[derive(Debug, Serialize)]
pub struct AncestralOutput {
  #[serde(skip)]
  pub gtr: Option<GTR>,
  pub model_name: GtrModelName,
  #[serde(skip)]
  pub mask: Vec<bool>,
  #[serde(skip)]
  pub emitted_nodes: Vec<GraphNodeKey>,
}
