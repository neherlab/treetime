use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
use crate::clock::reroot::reroot_undated;
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::optimize::dispatch::optimize_edge_lengths;
use crate::optimize::gather::InitialGuessEvidence;
use crate::optimize::indel::IndelPolicy;
use crate::optimize::iteration::apply_damping;
use crate::optimize::params::OptimizeLoopSettings;
use crate::optimize::params::{BranchOptMethod, InitialGuessMode};
use crate::optimize::partition_rates::{fix_equal_rates, normalize_partition_rates};
use crate::optimize::run_loop::{apply_initial_guess_mode, run_optimize_loop};
use crate::optimize::solver_outcome::ConvergenceSummary;
use crate::partition::alignment_sequences::assemble_sequence_mutations;
use crate::partition::create::Representation;
use crate::partition::marginal::reconstruction::ReconstructionRef;
use crate::partition::result::{
  AlignmentOutcome, AlignmentResult, AlignmentSequences, PartitionResult, alignment_results,
};
use crate::partition::set::{PartitionSet, build_partition_set};
use crate::partition::split::{AlignmentInfo, SequenceInput, SplitInput};
use crate::progress::{LogSink, StageSink};
use crate::seq::mutation::{MutationTrack, SequenceMutations};
use crate::{progress_info, progress_warn};
use eyre::Report;
use itertools::{Itertools, izip};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node_names::NodeNames;
use treetime_primitives::LogLh;
use treetime_utils::{make_internal_error, make_report};

const PRE_REROOT_DAMPING: f64 = 0.75;

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
  } = input;
  let PreparedPartitions {
    alignments,
    model_names,
    partitions,
  } = prepare_partitions(params, &graph, sequences, &mut branch_lengths, log).map_err(OperationError::classify)?;
  apply_initial_guess(params, &graph, &partitions, &mut branch_lengths, &names, log)
    .map_err(OperationError::classify)?;

  let mut convergence = ConvergenceSummary::default();
  let partitions = match &params.reroot_spec {
    Some(spec) => {
      progress_info!(log, "Rerooting before optimization: {spec:?}");
      stages.report("Rerooting", 0.2, "");
      let partitions = pre_reroot_optimize(
        &graph,
        partitions,
        params.loop_settings.opt_method,
        params.loop_settings.indel_policy,
        &mut branch_lengths,
        &names,
        &mut convergence,
      )
      .map_err(OperationError::classify)?;
      reroot_optimize(&mut graph, &mut names, spec, partitions, &mut branch_lengths, log)
        .map_err(OperationError::classify)?
    },
    None => partitions,
  };

  cancel.check().map_err(OperationError::classify)?;
  stages.report("Optimizing branch lengths", 0.3, "");
  let loop_result = run_optimize_loop(
    &mut graph,
    partitions,
    &params.loop_settings,
    branch_lengths,
    names,
    log,
  )
  .map_err(OperationError::classify)?;
  convergence.merge(loop_result.convergence);
  convergence.report(params.loop_settings.opt_method, log);
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

/// Builds the partitions of a tree, optimizes its branch lengths and partition rates without changing the
/// topology, and returns the partitions with the marginal messages on the optimized branch lengths
pub(crate) fn optimized_partitions(
  params: &OptimizeParams,
  graph: &mut Graph,
  names: NodeNames,
  sequences: SequenceInput,
  mut branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  log: &dyn LogSink,
) -> Result<(PartitionSet, BTreeMap<GraphEdgeKey, Option<f64>>, LogLh), Report> {
  let PreparedPartitions { partitions, .. } = prepare_partitions(params, graph, sequences, &mut branch_lengths, log)?;
  apply_initial_guess(params, graph, &partitions, &mut branch_lengths, &names, log)?;
  let result = run_optimize_loop(graph, partitions, &params.loop_settings, branch_lengths, names, log)?;
  let (partitions, log_lh) = result
    .partitions
    .marginal_update(graph, &branch_lengths_or_zero(&result.branch_lengths))?;
  Ok((partitions, result.branch_lengths, log_lh))
}

pub struct OptimizeInput {
  pub graph: Graph,
  pub names: NodeNames,
  pub sequences: SequenceInput,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

#[derive(Debug)]
pub struct OptimizeOutput {
  pub graph: Graph,
  pub alignments: Vec<AlignmentResult>,
  pub partitions: Vec<PartitionResult>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub names: NodeNames,
}

fn check_params(params: &OptimizeParams, sequences: &SequenceInput) -> Result<(), OperationError> {
  if !(0.0..1.0).contains(&params.loop_settings.damping) {
    return Err(OperationError::InvalidParams(make_report!(
      "damping must be in [0.0, 1.0), got {}",
      params.loop_settings.damping
    )));
  }

  if let Some(RerootSpec::Method(method)) = &params.reroot_spec
    && *method != RerootMethod::MinDev
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
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  log: &dyn LogSink,
) -> Result<PreparedPartitions, Report> {
  let normalize_inferred = sequences.has_inferred_model();
  let split = sequences.split()?;
  let model_names = split.model_names();
  let SplitInput { alignments, partitions } = split;
  let partitions = build_partition_set(
    Representation::resolve(params.dense),
    graph,
    partitions,
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
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
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
    .map(|edge| {
      let edge_key = edge.key();
      Ok((edge_key, reconstruction.edge_mutations(graph, edge_key, track)?))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  Ok(SequenceMutations {
    root_sequence: reconstruction.root_sequence(graph)?,
    edge_mutations,
  })
}

fn pre_reroot_optimize(
  graph: &Graph,
  partitions: PartitionSet,
  opt_method: BranchOptMethod,
  indel_policy: IndelPolicy,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &NodeNames,
  convergence: &mut ConvergenceSummary,
) -> Result<PartitionSet, Report> {
  let old_branch_lengths = branch_lengths.clone();

  let indels = indel_policy.edge_indels(graph, &partitions, branch_lengths)?;
  let unconverged = optimize_edge_lengths(
    graph,
    partitions.total_length(),
    &partitions,
    &indels,
    opt_method,
    branch_lengths,
  )?;
  convergence.record(&unconverged, graph, names);

  apply_damping(branch_lengths, &old_branch_lengths, PRE_REROOT_DAMPING, 0);
  Ok(
    partitions
      .marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?
      .0,
  )
}

fn reroot_optimize(
  graph: &mut Graph,
  names: &mut NodeNames,
  spec: &RerootSpec,
  partitions: PartitionSet,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  log: &dyn LogSink,
) -> Result<PartitionSet, Report> {
  if let RerootSpec::Method(method @ (RerootMethod::LeastSquares | RerootMethod::Oldest)) = spec {
    return make_internal_error!("optimize cannot reroot with a date-dependent method: {method:?}");
  }
  let reroot_result = reroot_undated(graph, names, branch_lengths, spec, log)?;

  let partitions = partitions.apply_reroot(&reroot_result)?;
  Ok(
    partitions
      .marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?
      .0,
  )
}
