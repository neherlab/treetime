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
use crate::optimize::run_loop::{apply_initial_guess_mode, normalize_partition_rates, run_optimize_loop};
use crate::optimize::solver_outcome::ConvergenceSummary;
use crate::partition::create::{Representation, build_marginal_partition};
use crate::partition::marginal::reconstruction::MarginalReconstruction;
use crate::partition::result::{PartitionResult, PartitionSequences, PartitionSource, SplitSource, split_source};
use crate::progress::{LogSink, StageSink};
use crate::reroot::orchestrate::{RerootTopologyParams, reroot_at_node, reroot_min_dev};
use crate::reroot::params::BrentParams;
use crate::reroot::variance::VarianceModel;
use crate::seq::alignment::SequenceSource;
use crate::seq::mutation::{MutationTrack, SequenceMutations};
use crate::{progress_info, progress_warn};
use eyre::Report;
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

  let mut branch_lengths = std::mem::take(&mut input.branch_lengths);
  let SplitSource {
    source,
    alphabet,
    nodes,
  } = split_source(1, input.source);

  let reconstruction = build_marginal_partition(
    Representation::resolve(params.dense),
    params.model,
    &input.graph,
    alphabet,
    nodes,
    &branch_lengths_or_zero(&branch_lengths),
    log,
  )
  .map_err(OperationError::classify)?;

  if matches!(reconstruction, MarginalReconstruction::Dense(_)) {
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

  let (mut reconstruction, _) = reconstruction
    .marginal_update(&input.graph, &branch_lengths_or_zero(&branch_lengths))
    .map_err(OperationError::classify)?;

  if params.model == GtrModelName::Infer {
    let length = reconstruction.sequence_length();
    normalize_partition_rates(&mut [(length, reconstruction.gtr_mut())], &mut branch_lengths)
      .map_err(OperationError::classify)?;
  }

  {
    let indels = params
      .loop_settings
      .indel_policy
      .edge_indels(&input.graph, &reconstruction, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&input.graph, &reconstruction).map_err(OperationError::classify)?;
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
    reconstruction = pre_reroot_optimize(
      &input.graph,
      reconstruction,
      params.loop_settings.opt_method,
      params.loop_settings.indel_policy,
      &mut branch_lengths,
      names,
      &mut convergence,
    )
    .map_err(OperationError::classify)?;
    reconstruction = reroot_optimize(&mut input.graph, spec, reconstruction, &mut branch_lengths, names)
      .map_err(OperationError::classify)?;
  }

  let loop_names = restrict_node_names(names, &input.graph);

  cancel.check().map_err(OperationError::classify)?;
  stages.report("Optimizing branch lengths", 0.3, "");
  let loop_result = run_optimize_loop(
    &mut input.graph,
    reconstruction,
    &params.loop_settings,
    branch_lengths,
    &loop_names,
  )
  .map_err(OperationError::classify)?;
  convergence.merge(loop_result.convergence);
  convergence.report(params.loop_settings.opt_method, log);
  let branch_lengths = loop_result.branch_lengths;

  progress_info!(log, "Re-running marginal to populate subs_ml after optimization loop");
  let (reconstruction, _) = loop_result
    .reconstruction
    .marginal_update(&input.graph, &branch_lengths_or_zero(&branch_lengths))
    .map_err(OperationError::classify)?;

  let result = optimize_result(
    &input.graph,
    source,
    &reconstruction,
    params.model,
    loop_result.indel_rate,
  )
  .map_err(OperationError::InferenceFailed)?;

  Ok(OptimizeOutput {
    graph: input.graph,
    results: vec![result],
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
  pub source: SequenceSource,
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
  reconstruction: &MarginalReconstruction,
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
    length: source.length,
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
  reconstruction: MarginalReconstruction,
  opt_method: BranchOptMethod,
  indel_policy: IndelPolicy,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  convergence: &mut ConvergenceSummary,
) -> Result<MarginalReconstruction, Report> {
  let old_branch_lengths = branch_lengths.clone();

  let indels = indel_policy.edge_indels(graph, &reconstruction, branch_lengths);
  let unconverged = optimize_edge_lengths(
    graph,
    reconstruction.sequence_length(),
    &reconstruction,
    &indels,
    opt_method,
    branch_lengths,
  )?;
  convergence.record(&unconverged, graph, names);

  apply_damping(branch_lengths, &old_branch_lengths, PRE_REROOT_DAMPING, 0);
  Ok(
    reconstruction
      .marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?
      .0,
  )
}

fn reroot_optimize(
  graph: &mut Graph,
  spec: &RerootSpec,
  reconstruction: MarginalReconstruction,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<MarginalReconstruction, Report> {
  let variance = VarianceModel::default();
  let topo = RerootTopologyParams::default();
  let opt_params = BrentParams::default();

  let reroot_result = match spec {
    RerootSpec::Method(RerootMethod::MinDev) => {
      reroot_min_dev(graph, &variance, &opt_params, topo, branch_lengths, names)?
    },
    RerootSpec::Tips(tips) => {
      let mrca = common_ancestor(graph, tips)?;
      reroot_at_node(graph, mrca, topo, branch_lengths, names)?
    },
    RerootSpec::Method(method) => {
      return make_internal_error!("optimize cannot reroot with a date-dependent method: {method:?}");
    },
  };

  let reconstruction = reconstruction.apply_reroot(&reroot_result)?;
  Ok(
    reconstruction
      .marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?
      .0,
  )
}
