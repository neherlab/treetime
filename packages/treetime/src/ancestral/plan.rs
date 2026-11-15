use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch::ancestral_reconstruction_fitch;
use crate::ancestral::params::AncestralParams;
use crate::ancestral::params::MethodAncestral;
use crate::ancestral::partition::AncestralPartition;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::gtr::gamma::NO_RATE_VARIATION;
use crate::gtr::get_gtr::{GtrModel, GtrModelName};
use crate::gtr::gtr::GTR;
use crate::optimize::partition_rates::scale_rates_to_unit_mean;
use crate::partition::create::{PartitionInput, Representation, build_marginal_partition, build_marginal_partitions};
use crate::partition::fitch::passes::create_inherited_fitch_partition;
use crate::partition::inferred_rate::InferredRate;
use crate::partition::marginal::reconstruction::MarginalReconstruction;
use crate::partition::marginal::sample::SampleMode;
use crate::partition::resolve::PartitionLayout;
use crate::partition::set::PartitionSet;
use crate::partition::split::{AlignmentInfo, ModelSettings, PartitionData, SplitInput, SplitPartition};
use crate::progress::{LogSink, StageSink};
use crate::seq::alignment::NodeSeqInput;
use crate::{progress_info, progress_warn};
use eyre::Report;
use itertools::{Itertools, izip};
use rand::Rng;
use std::sync::Arc;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::scope::{EdgeScopes, Inheritance};
use treetime_utils::{make_error, make_internal_report, make_report};

pub(crate) fn resolve_plan(params: &AncestralParams) -> Result<ReconstructionMethod, OperationError> {
  if params.sample_from_profile != SampleMode::Argmax && params.method != MethodAncestral::Marginal {
    return Err(OperationError::InvalidParams(make_report!(
      "--sample-from-profile={:?} requires --method-anc=marginal. Posterior sampling is only defined \
       for marginal reconstruction; {:?} has no posterior profile to sample.",
      params.sample_from_profile,
      params.method
    )));
  }

  match params.method {
    MethodAncestral::Parsimony => Ok(ReconstructionMethod::Fitch),
    MethodAncestral::Marginal => Ok(ReconstructionMethod::Marginal(Representation::resolve(params.dense))),
  }
}

#[expect(
  clippy::too_many_arguments,
  reason = "the partition inputs come from the command and the random generator, cancellation, stage and log sinks are run services of the caller; no type holds both groups"
)]
pub(crate) fn reconstruct_partitions(
  graph: &Graph,
  scopes: Option<&Arc<EdgeScopes>>,
  method: ReconstructionMethod,
  input: SplitInput,
  branch_lengths: &EdgeMap<f64>,
  options: &ReconstructionOptions,
  rng: &mut dyn Rng,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<Vec<ReconstructedPartition>, Report> {
  let SplitInput {
    alignments,
    partitions,
    models,
  } = input;
  let (heads, inputs): (Vec<_>, Vec<_>) = partitions
    .into_iter()
    .map(|SplitPartition { data, model }| {
      let PartitionData {
        layout,
        alphabet,
        nodes,
      } = data;
      ((layout, model), (alphabet, nodes))
    })
    .unzip();
  let inheritances = segment_inheritances(scopes, alignments.len(), &heads)?;
  let partitions = match method {
    ReconstructionMethod::Fitch => {
      let inputs = izip!(inheritances, inputs)
        .map(|(inheritance, (alphabet, nodes))| (inheritance, alphabet, nodes))
        .collect();
      reconstruct_fitch(graph, inputs, options.impute, cancel, stages, log)?
    },
    ReconstructionMethod::Marginal(_) if scopes.is_some() => {
      return make_error!(
        "Marginal reconstruction on a graph of several segment trees is not implemented yet; use --method-anc=parsimony"
      );
    },
    ReconstructionMethod::Marginal(representation) => {
      let reconstructions = build_alignment_partitions(
        graph,
        representation,
        &alignments,
        &heads,
        inputs,
        models,
        branch_lengths,
        cancel,
        stages,
        log,
      )?;
      let reconstructions = estimate_gamma_shapes(reconstructions, graph, branch_lengths, log)?;
      reconstructions
        .into_iter()
        .map(|reconstruction| {
          checkpoint(cancel, stages, "Marginal reconstruction", 0.4)?;
          reconstruct_marginal(graph, reconstruction, branch_lengths, options, rng, cancel, stages, log)
        })
        .collect::<Result<Vec<_>, Report>>()?
    },
  };
  Ok(
    izip!(heads, partitions)
      .map(|((layout, model), partition)| ReconstructedPartition {
        layout,
        model_name: model.name(),
        partition,
      })
      .collect(),
  )
}

fn estimate_gamma_shapes(
  reconstructions: Vec<MarginalReconstruction>,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  log: &dyn LogSink,
) -> Result<Vec<MarginalReconstruction>, Report> {
  if reconstructions
    .iter()
    .all(|reconstruction| reconstruction.model().categories() == 1)
  {
    return Ok(reconstructions);
  }
  Ok(
    PartitionSet::from_reconstructions(reconstructions)?
      .estimate_gamma_shapes(graph, branch_lengths, log)?
      .into_reconstructions(),
  )
}

fn segment_inheritances(
  scopes: Option<&Arc<EdgeScopes>>,
  alignment_count: usize,
  heads: &[(PartitionLayout, GtrModel)],
) -> Result<Vec<Inheritance>, Report> {
  let Some(scopes) = scopes else {
    return Ok(vec![Inheritance::ALL; heads.len()]);
  };
  if scopes.segment_count() != alignment_count {
    return make_error!(
      "The graph has {} segment trees, but the run has {} alignments: give one alignment per segment tree, in the same order",
      scopes.segment_count(),
      alignment_count
    );
  }
  heads
    .iter()
    .map(|(layout, _)| Inheritance::segment(Arc::clone(scopes), layout.alignment - 1))
    .collect()
}

pub(crate) struct ReconstructedPartition {
  pub(crate) layout: PartitionLayout,
  pub(crate) model_name: GtrModelName,
  pub(crate) partition: AncestralPartition,
}

#[expect(
  clippy::too_many_arguments,
  reason = "the partition inputs come from the command and the cancellation, stage and log sinks are run services of the caller; no type holds both groups"
)]
fn build_alignment_partitions(
  graph: &Graph,
  representation: Representation,
  alignments: &[AlignmentInfo],
  heads: &[(PartitionLayout, GtrModel)],
  inputs: Vec<(Alphabet, NodeMap<NodeSeqInput>)>,
  models: ModelSettings,
  branch_lengths: &EdgeMap<f64>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<Vec<MarginalReconstruction>, Report> {
  let counts = heads.iter().counts_by(|(layout, _)| layout.alignment);
  checkpoint(cancel, stages, "Inferring GTR model", 0.2)?;
  if heads.len() > 1 {
    for (layout, _) in heads {
      progress_info!(
        log,
        "{}: alignment '{}'",
        capitalized_label(layout),
        alignments[layout.alignment - 1].path.display()
      );
    }
  }
  let inputs = izip!(heads, inputs)
    .map(|((layout, model), (alphabet, nodes))| PartitionInput {
      model,
      rate: InferredRate::for_partition_count(counts[&layout.alignment]),
      alphabet,
      nodes,
    })
    .collect_vec();
  let mut reconstructions = izip!(
    heads,
    build_marginal_partitions(representation, inputs, models, graph, branch_lengths, log)?
  )
  .map(|((layout, _), reconstruction)| {
    reconstruction
      .with_layout(layout.clone())
      .with_total_length(alignments[layout.alignment - 1].length)
  })
  .collect_vec();
  for alignment in alignments {
    let mut inferred = izip!(heads, &mut reconstructions)
      .filter(|((layout, model), _)| layout.alignment == alignment.number && model.is_infer())
      .map(|((layout, _), reconstruction)| (layout.columns.len(), reconstruction_gtr(reconstruction)))
      .collect_vec();
    if counts[&alignment.number] > 1 {
      scale_rates_to_unit_mean(&mut inferred)?;
    }
  }
  Ok(reconstructions)
}

fn reconstruction_gtr(reconstruction: &mut MarginalReconstruction) -> &mut GTR {
  match reconstruction {
    MarginalReconstruction::Dense(dense) => dense.model.gtr_mut(),
    MarginalReconstruction::Sparse(sparse) => sparse.model.gtr_mut(),
  }
}

fn capitalized_label(layout: &PartitionLayout) -> String {
  match &layout.name {
    Some(name) => format!("Partition {} '{name}'", layout.number),
    None => format!("Partition {}", layout.number),
  }
}

#[expect(
  clippy::too_many_arguments,
  reason = "the partition inputs come from the command and the random generator, cancellation, stage and log sinks are run services of the caller; no type holds both groups"
)]
pub(crate) fn reconstruct_partition(
  graph: &Graph,
  method: ReconstructionMethod,
  model: &GtrModel,
  rate: InferredRate,
  alphabet: Alphabet,
  node_inputs: NodeMap<NodeSeqInput>,
  branch_lengths: &EdgeMap<f64>,
  options: &ReconstructionOptions,
  rng: &mut dyn Rng,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<AncestralPartition, Report> {
  match method {
    ReconstructionMethod::Fitch => reconstruct_fitch(
      graph,
      vec![(Inheritance::ALL, alphabet, node_inputs)],
      options.impute,
      cancel,
      stages,
      log,
    )?
    .pop()
    .ok_or_else(|| make_internal_report!("Fitch reconstruction lost its partition")),
    ReconstructionMethod::Marginal(representation) => {
      checkpoint(cancel, stages, "Inferring GTR model", 0.2)?;
      let reconstruction = build_marginal_partition(
        representation,
        model,
        rate,
        NO_RATE_VARIATION,
        graph,
        alphabet,
        node_inputs,
        branch_lengths,
        log,
      )?;
      checkpoint(cancel, stages, "Marginal reconstruction", 0.4)?;
      reconstruct_marginal(graph, reconstruction, branch_lengths, options, rng, cancel, stages, log)
    },
  }
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum ReconstructionMethod {
  Fitch,
  Marginal(Representation),
}

fn reconstruct_fitch(
  graph: &Graph,
  inputs: Vec<(Inheritance, Alphabet, NodeMap<NodeSeqInput>)>,
  impute: bool,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<Vec<AncestralPartition>, Report> {
  checkpoint(cancel, stages, "Fitch parsimony", 0.3)?;
  let mut partitions = inputs
    .into_iter()
    .map(|(inheritance, alphabet, node_inputs)| {
      create_inherited_fitch_partition(graph, inheritance, alphabet, node_inputs)
    })
    .collect::<Result<Vec<_>, Report>>()?;

  if impute {
    progress_warn!(
      log,
      "--impute-missing-data has no effect with --method-anc=parsimony: Fitch parsimony produces no \
       posterior profile to impute missing tip states from. Leaf states are emitted as observed."
    );
  }

  ancestral_reconstruction_fitch(graph, &mut partitions)?;
  Ok(partitions.into_iter().map(AncestralPartition::Fitch).collect())
}

#[expect(
  clippy::too_many_arguments,
  reason = "the partition inputs come from the command and the random generator, cancellation, stage and log sinks are run services of the caller; no type holds both groups"
)]
pub(crate) fn reconstruct_marginal(
  graph: &Graph,
  reconstruction: MarginalReconstruction,
  branch_lengths: &EdgeMap<f64>,
  options: &ReconstructionOptions,
  rng: &mut dyn Rng,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<AncestralPartition, Report> {
  let reconstruction = fit_reconstruction(reconstruction, graph, branch_lengths, options.impute)?;
  checkpoint(cancel, stages, "Reconstructing sequences", 0.6)?;
  let sampled = reconstruction.sample_sequences(graph, options.sample_mode, rng)?;
  Ok(AncestralPartition::Marginal {
    reconstruction,
    sampled,
    impute: options.impute,
  })
}

pub(crate) struct ReconstructionOptions {
  pub(crate) impute: bool,
  pub(crate) sample_mode: SampleMode,
}

impl ReconstructionOptions {
  pub(crate) const fn new(impute: bool, sample_mode: SampleMode) -> Self {
    Self { impute, sample_mode }
  }
}

fn fit_reconstruction(
  reconstruction: MarginalReconstruction,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  impute: bool,
) -> Result<MarginalReconstruction, Report> {
  Ok(match reconstruction {
    MarginalReconstruction::Sparse(sparse) => {
      MarginalReconstruction::Sparse(sparse.with_imputation(impute).marginal_update(graph, branch_lengths)?.0)
    },
    MarginalReconstruction::Dense(dense) => {
      MarginalReconstruction::Dense(dense.marginal_update(graph, branch_lengths)?.0)
    },
  })
}

fn checkpoint(cancel: &dyn Cancel, stages: &dyn StageSink, stage: &str, fraction: f64) -> Result<(), Report> {
  cancel.check()?;
  stages.report(stage, fraction, "");
  Ok(())
}
