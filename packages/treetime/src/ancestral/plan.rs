use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch::ancestral_reconstruction_fitch;
use crate::ancestral::params::AncestralParams;
use crate::ancestral::params::MethodAncestral;
use crate::ancestral::partition::AncestralPartition;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::gtr::refinement::refine_gtr_model;
use crate::partition::create::{Representation, build_marginal_partition};
use crate::partition::fitch::passes::create_fitch_partition;
use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction, SparseReconstruction};
use crate::partition::marginal::sample::SampleMode;
use crate::partition::marginal::shared::update::{MarginalPasses, MarginalUpdate};
use crate::progress::{LogSink, StageSink};
use crate::progress_warn;
use crate::seq::alignment::NodeSeqInput;
use eyre::Report;
use rand::RngCore;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::{make_internal_report, make_report};

#[derive(Clone, Copy, Debug)]
pub(crate) enum ReconstructionPlan {
  Fitch,
  Marginal {
    representation: Representation,
    model: GtrModelName,
    gtr_refinement: Option<usize>,
  },
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

pub(crate) fn resolve_plan(params: &AncestralParams) -> Result<ReconstructionPlan, OperationError> {
  if params.site_specific_gtr {
    return Err(OperationError::InvalidParams(make_report!(
      "--site-specific-gtr is not implemented"
    )));
  }

  if params.sample_from_profile != SampleMode::Argmax && params.method != MethodAncestral::Marginal {
    return Err(OperationError::InvalidParams(make_report!(
      "--sample-from-profile={:?} requires --method-anc=marginal. Posterior sampling is only defined \
       for marginal reconstruction; {:?} has no posterior profile to sample.",
      params.sample_from_profile,
      params.method
    )));
  }

  match params.method {
    MethodAncestral::Parsimony => Ok(ReconstructionPlan::Fitch),
    MethodAncestral::Marginal => Ok(ReconstructionPlan::Marginal {
      representation: Representation::resolve(params.dense),
      model: params.model,
      gtr_refinement: (params.gtr_iterations > 0 && params.model == GtrModelName::Infer)
        .then_some(params.gtr_iterations),
    }),
  }
}

pub(crate) fn reconstruct_partition(
  graph: &Graph,
  plan: &ReconstructionPlan,
  index: usize,
  alphabet: Alphabet,
  node_inputs: BTreeMap<GraphNodeKey, NodeSeqInput>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  options: &ReconstructionOptions,
  rng: &mut dyn RngCore,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<AncestralPartition, Report> {
  match *plan {
    ReconstructionPlan::Fitch => {
      reconstruct_fitch(graph, index, alphabet, node_inputs, options.impute, cancel, stages, log)
    },
    ReconstructionPlan::Marginal {
      representation,
      model,
      gtr_refinement,
    } => {
      checkpoint(cancel, stages, "Inferring GTR model", 0.2)?;
      let reconstruction = build_marginal_partition(
        representation,
        model,
        graph,
        index,
        alphabet,
        node_inputs,
        branch_lengths,
        log,
      )?;
      checkpoint(cancel, stages, "Marginal reconstruction", 0.4)?;
      reconstruct_marginal(
        graph,
        reconstruction,
        gtr_refinement,
        branch_lengths,
        options,
        rng,
        cancel,
        stages,
        log,
      )
    },
  }
}

fn reconstruct_fitch(
  graph: &Graph,
  index: usize,
  alphabet: Alphabet,
  node_inputs: BTreeMap<GraphNodeKey, NodeSeqInput>,
  impute: bool,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<AncestralPartition, Report> {
  checkpoint(cancel, stages, "Fitch parsimony", 0.3)?;
  let mut partitions = vec![create_fitch_partition(graph, index, alphabet, node_inputs)?];

  if impute {
    progress_warn!(
      log,
      "--impute-missing-data has no effect with --method-anc=parsimony: Fitch parsimony produces no \
       posterior profile to impute missing tip states from. Leaf states are emitted as observed."
    );
  }

  ancestral_reconstruction_fitch(graph, &mut partitions)?;
  let partition = partitions
    .pop()
    .ok_or_else(|| make_internal_report!("Fitch reconstruction lost its partition"))?;
  Ok(AncestralPartition::Fitch(partition))
}

fn reconstruct_marginal(
  graph: &Graph,
  reconstruction: MarginalReconstruction,
  gtr_refinement: Option<usize>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  options: &ReconstructionOptions,
  rng: &mut dyn RngCore,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<AncestralPartition, Report> {
  let reconstruction = match gtr_refinement {
    None => reconstruction.marginal_update(graph, branch_lengths)?.0,
    Some(iterations) => refine_gtr(reconstruction, iterations, graph, branch_lengths, log)?,
  };
  checkpoint(cancel, stages, "Reconstructing sequences", 0.6)?;
  let sampled = reconstruction.sample_sequences(graph, options.sample_mode, rng)?;
  Ok(AncestralPartition::Marginal {
    reconstruction,
    sampled,
    impute: options.impute,
  })
}

fn refine_gtr(
  reconstruction: MarginalReconstruction,
  iterations: usize,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<MarginalReconstruction, Report> {
  Ok(match reconstruction {
    MarginalReconstruction::Sparse(SparseReconstruction {
      partition,
      gtr,
      node_states,
      ..
    }) => {
      let (gtr, MarginalUpdate { node_states, edges, .. }) =
        refined_update(&partition, gtr, &node_states, iterations, graph, branch_lengths, log)?;
      MarginalReconstruction::Sparse(SparseReconstruction {
        partition,
        gtr,
        node_states,
        edges,
      })
    },
    MarginalReconstruction::Dense(DenseReconstruction { partition, gtr, .. }) => {
      let (gtr, MarginalUpdate { node_states, edges, .. }) =
        refined_update(&partition, gtr, &(), iterations, graph, branch_lengths, log)?;
      MarginalReconstruction::Dense(DenseReconstruction {
        partition,
        gtr,
        node_states,
        edges,
      })
    },
  })
}

fn refined_update<P: MarginalPasses>(
  partition: &P,
  gtr: GTR,
  input: &P::BackwardInput,
  iterations: usize,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  log: &dyn LogSink,
) -> Result<(GTR, MarginalUpdate<P::Node, P::Backward, P::Forward, P::Estimate>), Report> {
  let update = partition.marginal_update(&gtr, graph, branch_lengths, input)?;
  refine_gtr_model(partition, gtr, update, iterations, 1.0, graph, branch_lengths, log)
}

fn checkpoint(cancel: &dyn Cancel, stages: &dyn StageSink, stage: &str, fraction: f64) -> Result<(), Report> {
  cancel.check()?;
  stages.report(stage, fraction, "");
  Ok(())
}
