use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
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
use crate::partition::create::Representation;
use crate::partition::marginal::reconstruction::ReconstructionRef;
use crate::partition::result::{PartitionResult, PartitionSequences, PartitionSource};
use crate::partition::set::{PartitionSet, build_partition_set};
use crate::progress::{LogSink, StageSink};
use crate::reroot::orchestrate::{RerootTopologyParams, reroot_at_node, reroot_min_dev};
use crate::seq::alignment::SequenceSource;
use crate::seq::mutation::{MutationTrack, SequenceMutations};
use crate::{progress_info, progress_warn};
use eyre::Report;
use itertools::izip;
use std::collections::BTreeMap;
use treetime_graph::assign_node_names::restrict_node_names;
use treetime_graph::common_ancestor::common_ancestor;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::{make_internal_error, make_report};

const PRE_REROOT_DAMPING: f64 = 0.75;

pub fn run(
  params: &OptimizeParams,
  mut input: OptimizeInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<OptimizeOutput, OperationError> {
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

  if input.sources.is_empty() {
    return Err(OperationError::InvalidInput(make_report!(
      "Branch length optimization requires at least one alignment"
    )));
  }

  let mut branch_lengths = std::mem::take(&mut input.branch_lengths);
  let (sources, partitions) = build_partition_set(
    Representation::resolve(params.dense),
    params.model,
    &input.graph,
    std::mem::take(&mut input.sources),
    &branch_lengths_or_zero(&branch_lengths),
    log,
  )
  .map_err(OperationError::classify)?;

  if matches!(partitions, PartitionSet::Dense(_)) {
    if !params.loop_settings.topology_ops.merge_siblings {
      progress_warn!(
        log,
        "--no-merge-siblings has no effect in dense mode: sibling merging requires the sparse sequence representation"
      );
    }
    if !params.loop_settings.topology_ops.flip_parent_child {
      progress_warn!(
        log,
        "--no-flip-parent-child has no effect in dense mode: the reversion hoist requires the sparse sequence representation"
      );
    }
  }

  let (mut partitions, _) = partitions
    .marginal_update(&input.graph, &branch_lengths_or_zero(&branch_lengths))
    .map_err(OperationError::classify)?;

  if params.model == GtrModelName::Infer {
    normalize_partition_rates(&mut partitions, &mut branch_lengths, None).map_err(OperationError::classify)?;
  }
  let mut partitions = fix_equal_rates(
    &input.graph,
    partitions,
    params.loop_settings.partition_rates,
    &branch_lengths,
  )
  .map_err(OperationError::classify)?;

  {
    let indels = params
      .loop_settings
      .indel_policy
      .edge_indels(&input.graph, &partitions, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&input.graph, &partitions).map_err(OperationError::classify)?;
    apply_initial_guess_mode(
      &input.graph,
      &indels,
      &evidence,
      params.initial_guess,
      &mut branch_lengths,
      names,
      log,
    )
    .map_err(OperationError::classify)?;
  }

  let mut convergence = ConvergenceSummary::default();
  if let Some(spec) = &params.reroot_spec {
    progress_info!(log, "Rerooting before optimization: {spec:?}");
    stages.report("Rerooting", 0.2, "");
    partitions = pre_reroot_optimize(
      &input.graph,
      partitions,
      params.loop_settings.opt_method,
      params.loop_settings.indel_policy,
      &mut branch_lengths,
      names,
      &mut convergence,
    )
    .map_err(OperationError::classify)?;
    partitions = reroot_optimize(&mut input.graph, spec, partitions, &mut branch_lengths, names, log)
      .map_err(OperationError::classify)?;
  }

  let loop_names = restrict_node_names(names, &input.graph);

  cancel.check().map_err(OperationError::classify)?;
  stages.report("Optimizing branch lengths", 0.3, "");
  let loop_result = run_optimize_loop(
    &mut input.graph,
    partitions,
    &params.loop_settings,
    branch_lengths,
    &loop_names,
    log,
  )
  .map_err(OperationError::classify)?;
  convergence.merge(loop_result.convergence);
  convergence.report(params.loop_settings.opt_method, log);
  let branch_lengths = loop_result.branch_lengths;

  progress_info!(log, "Re-running marginal to populate subs_ml after optimization loop");
  let (partitions, _) = loop_result
    .partitions
    .marginal_update(&input.graph, &branch_lengths_or_zero(&branch_lengths))
    .map_err(OperationError::classify)?;

  let results = izip!(sources, partitions.iter(), loop_result.indel_rates)
    .map(|(source, partition, indel_rate)| optimize_result(&input.graph, source, partition, params.model, indel_rate))
    .collect::<Result<Vec<_>, Report>>()
    .map_err(OperationError::InferenceFailed)?;

  Ok(OptimizeOutput {
    graph: input.graph,
    results,
    branch_lengths,
    names: loop_result.names,
  })
}

pub struct OptimizeParams {
  pub model: GtrModelName,
  pub dense: Option<bool>,
  pub initial_guess: InitialGuessMode,
  pub reroot_spec: Option<RerootSpec>,
  pub loop_settings: OptimizeLoopSettings,
}

pub struct OptimizeInput {
  pub graph: Graph,
  pub sources: Vec<SequenceSource>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

#[derive(Debug)]
pub struct OptimizeOutput {
  pub graph: Graph,
  pub results: Vec<PartitionResult>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
}

fn optimize_result(
  graph: &Graph,
  source: PartitionSource,
  reconstruction: ReconstructionRef<'_>,
  model_name: GtrModelName,
  indel_rate: Option<f64>,
) -> Result<PartitionResult, Report> {
  let edge_mutations = graph
    .get_edges()
    .map(|edge| {
      let edge_key = edge.key();
      Ok((
        edge_key,
        reconstruction.edge_mutations(graph, edge_key, &MutationTrack::Nucleotide)?,
      ))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  let mutations = SequenceMutations {
    root_sequence: reconstruction.root_sequence(graph)?,
    edge_mutations,
  };
  let alphabet = reconstruction.alphabet().clone();
  Ok(PartitionResult {
    number: source.number,
    path: source.path,
    length: reconstruction.length(),
    mask: source.mask,
    model_name,
    gtr: Some(reconstruction.gtr().clone()),
    indel_rate,
    log_lh: Some(reconstruction.graph_log_lh(graph)?),
    sequences: Some(PartitionSequences::new(mutations, &alphabet)?),
    alphabet,
  })
}

fn pre_reroot_optimize(
  graph: &Graph,
  partitions: PartitionSet,
  opt_method: BranchOptMethod,
  indel_policy: IndelPolicy,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  convergence: &mut ConvergenceSummary,
) -> Result<PartitionSet, Report> {
  let old_branch_lengths = branch_lengths.clone();

  let indels = indel_policy.edge_indels(graph, &partitions, branch_lengths);
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
  spec: &RerootSpec,
  partitions: PartitionSet,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<PartitionSet, Report> {
  let topo = RerootTopologyParams::default();

  let reroot_result = match spec {
    RerootSpec::Method(RerootMethod::MinDev) => reroot_min_dev(graph, topo, branch_lengths, names, log)?,
    RerootSpec::Tips(tips) => {
      let mrca = common_ancestor(graph, tips)?;
      reroot_at_node(graph, mrca, topo, branch_lengths, names)?
    },
    RerootSpec::Method(method) => {
      return make_internal_error!("optimize cannot reroot with a date-dependent method: {method:?}");
    },
  };

  let partitions = partitions.apply_reroot(&reroot_result)?;
  Ok(
    partitions
      .marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?
      .0,
  )
}
