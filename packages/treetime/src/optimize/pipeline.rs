use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::optimize::gather::InitialGuessEvidence;
use crate::optimize::params::InitialGuessMode;
use crate::optimize::params::{OptimizeLoopSettings, TopologyOps};
use crate::optimize::partition_rates::{fix_equal_rates, normalize_partition_rates};
use crate::optimize::run_loop::{apply_initial_guess_mode, run_optimize_loop};
use crate::partition::alignment_sequences::assemble_sequence_mutations;
use crate::partition::create::{Representation, segment_inheritances};
use crate::partition::marginal::reconstruction::ReconstructionRef;
use crate::partition::result::{
  AlignmentOutcome, AlignmentResult, AlignmentSequences, PartitionResult, alignment_results,
};
use crate::partition::set::{PartitionSet, build_inherited_partition_set};
use crate::partition::split::{AlignmentInfo, SequenceInput, SplitInput};
use crate::progress::{LogSink, StageSink};
use crate::reroot::date_free::{LengthsAdvice, reroot_date_free};
use crate::reroot::params::RerootSpec;
use crate::seq::mutation::{MutationTrack, SequenceMutations};
use crate::{progress_info, progress_warn};
use eyre::Report;
use itertools::{Itertools, izip};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::node_names::NodeNames;
use treetime_graph::scope::EdgeScopes;
use treetime_utils::make_report;

pub fn run(
  params: &OptimizeParams,
  input: OptimizeInput,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<OptimizeOutput, OperationError> {
  check_params(params, &input.sequences)?;
  let OptimizeInput {
    mut graph,
    mut names,
    sequences,
    mut branch_lengths,
    scopes,
  } = input;
  let loop_settings = segment_loop_settings(params, scopes.as_ref(), log)?;
  let PreparedPartitions {
    alignments,
    model_names,
    partitions,
  } = prepare_partitions(params, &graph, sequences, scopes.as_ref(), &mut branch_lengths, log)
    .map_err(OperationError::classify)?;
  apply_initial_guess(params, &graph, &partitions, &mut branch_lengths, &names, log)
    .map_err(OperationError::classify)?;

  let partitions = match &params.reroot_spec {
    Some(spec) => {
      progress_info!(log, "Rerooting before optimization: {spec:?}");
      stages.report("Rerooting", 0.2, "");
      reroot_optimize(&mut graph, &mut names, spec, partitions, &mut branch_lengths, log)
        .map_err(OperationError::classify)?
    },
    None => partitions,
  };

  cancel.check().map_err(OperationError::classify)?;
  stages.report("Optimizing branch lengths", 0.3, "");
  let loop_result = run_optimize_loop(&mut graph, partitions, &loop_settings, branch_lengths, names, log)
    .map_err(OperationError::classify)?;
  loop_result.convergence.report(loop_settings.opt_method, log);
  let branch_lengths = loop_result.branch_lengths;

  progress_info!(log, "Re-running marginal to populate subs_ml after optimization loop");
  let (partitions, _) = loop_result
    .partitions
    .marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))
    .map_err(OperationError::classify)?;

  let (alignments, partitions) =
    optimize_results(&graph, alignments, &partitions, &model_names, &loop_result.indel_rates)
      .map_err(OperationError::InferenceFailed)?;

  Ok(OptimizeOutput {
    graph,
    alignments,
    partitions,
    branch_lengths,
    names: loop_result.names,
  })
}

pub struct OptimizeInput {
  pub graph: Graph,
  pub names: NodeNames,
  pub sequences: SequenceInput,
  pub branch_lengths: EdgeMap<Option<f64>>,
  pub scopes: Option<Arc<EdgeScopes>>,
}

fn segment_loop_settings(
  params: &OptimizeParams,
  scopes: Option<&Arc<EdgeScopes>>,
  log: &dyn LogSink,
) -> Result<OptimizeLoopSettings, OperationError> {
  let mut settings = params.loop_settings;
  if scopes.is_none() {
    return Ok(settings);
  }
  if params.reroot_spec.is_some() {
    return Err(OperationError::InvalidParams(make_report!(
      "optimize cannot reroot a graph of several segment trees"
    )));
  }
  if settings.topology_ops != TopologyOps::NONE {
    progress_warn!(
      log,
      "Topology edits (zero-length collapse, merging of shared mutations, reversion hoisting) are off on a graph of several segment trees"
    );
    settings.topology_ops = TopologyOps::NONE;
  }
  Ok(settings)
}

#[derive(Debug)]
pub struct OptimizeOutput {
  pub graph: Graph,
  pub alignments: Vec<AlignmentResult>,
  pub partitions: Vec<PartitionResult>,
  pub branch_lengths: EdgeMap<Option<f64>>,
  pub names: NodeNames,
}

fn check_params(params: &OptimizeParams, sequences: &SequenceInput) -> Result<(), OperationError> {
  if !(0.0..1.0).contains(&params.loop_settings.damping) {
    return Err(OperationError::InvalidParams(make_report!(
      "damping must be in [0.0, 1.0), got {}",
      params.loop_settings.damping
    )));
  }

  if let Some(spec @ RerootSpec::Method(method)) = &params.reroot_spec
    && spec.uses_dates()
  {
    return Err(OperationError::InvalidParams(make_report!(
      "optimize cannot reroot with a date-dependent method: {method:?}"
    )));
  }

  if sequences.is_empty() {
    return Err(OperationError::InvalidInput(make_report!(
      "Branch length optimization requires at least one alignment"
    )));
  }
  Ok(())
}

fn prepare_partitions(
  params: &OptimizeParams,
  graph: &Graph,
  sequences: SequenceInput,
  scopes: Option<&Arc<EdgeScopes>>,
  branch_lengths: &mut EdgeMap<Option<f64>>,
  log: &dyn LogSink,
) -> Result<PreparedPartitions, Report> {
  let normalize_inferred = sequences.has_inferred_model();
  let split = sequences.split()?;
  let model_names = split.model_names();
  let SplitInput { alignments, partitions } = split;
  let inheritances = segment_inheritances(
    scopes,
    alignments.len(),
    partitions.iter().map(|partition| &partition.data.layout),
  )?;
  let partitions = build_inherited_partition_set(
    Representation::resolve(params.dense),
    graph,
    partitions,
    inheritances,
    &branch_lengths_or_zero(branch_lengths),
    log,
  )?;

  if matches!(partitions, PartitionSet::Dense(_)) {
    warn_sparse_only_topology_ops(params, log);
  }

  let (mut partitions, _) = partitions.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;

  if normalize_inferred {
    normalize_partition_rates(&mut partitions, branch_lengths, None)?;
  }
  let partitions = fix_equal_rates(graph, partitions, params.loop_settings.partition_rates, branch_lengths)?;
  Ok(PreparedPartitions {
    alignments,
    model_names,
    partitions,
  })
}

struct PreparedPartitions {
  alignments: Vec<AlignmentInfo>,
  model_names: Vec<GtrModelName>,
  partitions: PartitionSet,
}

fn warn_sparse_only_topology_ops(params: &OptimizeParams, log: &dyn LogSink) {
  if !params.loop_settings.topology_ops.merge_siblings {
    progress_warn!(
      log,
      "--no-merge-shared-mutations has no effect with --representation dense: sibling merging requires the sparse sequence representation"
    );
  }
  if !params.loop_settings.topology_ops.flip_parent_child {
    progress_warn!(
      log,
      "--no-hoist-reversions has no effect with --representation dense: the reversion hoist requires the sparse sequence representation"
    );
  }
}

fn apply_initial_guess(
  params: &OptimizeParams,
  graph: &Graph,
  partitions: &PartitionSet,
  branch_lengths: &mut EdgeMap<Option<f64>>,
  names: &NodeNames,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let indels = params
    .loop_settings
    .indel_policy
    .edge_indels(graph, partitions, branch_lengths)?;
  let evidence = InitialGuessEvidence::gather(graph, partitions)?;
  apply_initial_guess_mode(
    graph,
    &indels,
    &evidence,
    params.initial_guess,
    branch_lengths,
    names,
    log,
  )
}

pub struct OptimizeParams {
  pub dense: Option<bool>,
  pub initial_guess: InitialGuessMode,
  pub reroot_spec: Option<RerootSpec>,
  pub loop_settings: OptimizeLoopSettings,
}

fn optimize_results(
  graph: &Graph,
  alignments: Vec<AlignmentInfo>,
  partitions: &PartitionSet,
  model_names: &[GtrModelName],
  indel_rates: &BTreeMap<usize, f64>,
) -> Result<(Vec<AlignmentResult>, Vec<PartitionResult>), Report> {
  let reconstructions = partitions.iter().collect_vec();
  let track = MutationTrack::Nucleotide;
  let outcomes = partitions
    .alignment_groups()?
    .iter()
    .map(|group| {
      let members = group.members(&reconstructions);
      let mutations = members
        .iter()
        .map(|reconstruction| partition_mutations(graph, **reconstruction, &track))
        .try_collect()?;
      let mutations = assemble_sequence_mutations(
        &group.layout,
        mutations,
        |member, edge_key| members[member].edge_indels(edge_key),
        &track,
      )?;
      Ok(AlignmentOutcome {
        alignment: group.alignment(),
        indel_rate: indel_rates.get(&group.alignment()).copied(),
        sequences: Some(AlignmentSequences::new(mutations, members[0].alphabet())?),
      })
    })
    .collect::<Result<Vec<_>, Report>>()?;
  let results = izip!(partitions.iter(), partitions.layouts(), model_names)
    .map(|(partition, layout, &model_name)| {
      Ok(PartitionResult::new(
        &layout,
        model_name,
        Some(partition.gtr().clone()),
        Some(partition.graph_log_lh(graph)?),
      ))
    })
    .collect::<Result<Vec<_>, Report>>()?;
  Ok((alignment_results(alignments, outcomes)?, results))
}

fn partition_mutations(
  graph: &Graph,
  reconstruction: ReconstructionRef<'_>,
  track: &MutationTrack,
) -> Result<SequenceMutations, Report> {
  let edge_mutations = graph
    .get_edges()
    .filter(|edge| reconstruction.carries(edge.key()))
    .map(|edge| {
      let edge_key = edge.key();
      Ok((edge_key, reconstruction.edge_mutations(graph, edge_key, track)?))
    })
    .collect::<Result<EdgeMap<_>, Report>>()?;
  Ok(SequenceMutations {
    root_sequence: reconstruction.root_sequence(graph)?,
    edge_mutations,
  })
}

fn reroot_optimize(
  graph: &mut Graph,
  names: &mut NodeNames,
  spec: &RerootSpec,
  partitions: PartitionSet,
  branch_lengths: &mut EdgeMap<Option<f64>>,
  log: &dyn LogSink,
) -> Result<PartitionSet, Report> {
  let reroot_result = reroot_date_free(
    graph,
    branch_lengths,
    spec,
    &|_| false,
    LengthsAdvice::InitialBranchLengths,
    names,
    log,
  )?;

  let partitions = partitions.apply_reroot(&reroot_result)?;
  Ok(
    partitions
      .marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?
      .0,
  )
}
