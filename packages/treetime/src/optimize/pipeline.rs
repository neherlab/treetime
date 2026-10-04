use crate::alphabet::alphabet::Alphabet;
use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
use crate::error::{OperationError, input_error};
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::optimize::dispatch::{run_optimize_mixed, run_optimize_mixed_inner};
use crate::optimize::gather::{
  gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts,
};
use crate::optimize::iteration::apply_damping;
use crate::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};
use crate::optimize::run_loop::{apply_initial_guess_mode, normalize_partition_rates, run_optimize_loop};
use crate::partition::create::{Representation, build_marginal_partition};
use crate::partition::marginal::reconstruction::MarginalReconstruction;
use crate::progress::{LogSink, StageSink};
use crate::reroot::orchestrate::{RerootTopologyParams, reroot_at_node, reroot_min_dev};
use crate::reroot::params::BrentParams;
use crate::reroot::variance::VarianceModel;
use crate::seq::alignment::node_seq_inputs;
use crate::{progress_info, progress_warn};
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::assign_node_names::restrict_node_names;
use treetime_graph::common_ancestor::common_ancestor;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::AlignmentRecord;
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
  if !(0.0..1.0).contains(&params.damping) {
    return Err(OperationError::InvalidParams(make_report!(
      "damping must be in [0.0, 1.0), got {}",
      params.damping
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
  let sequences = std::mem::take(&mut input.sequences);
  let node_inputs = node_seq_inputs(&input.graph, names, sequences);

  let reconstruction = build_marginal_partition(
    Representation::resolve(params.dense),
    params.model,
    &input.graph,
    0,
    input.alphabet,
    &node_inputs,
    &branch_lengths_or_zero(&branch_lengths),
    log,
  )
  .map_err(OperationError::classify)?;

  if matches!(reconstruction, MarginalReconstruction::Dense(_)) {
    if !params.topology_ops.merge_siblings {
      progress_warn!(
        log,
        "--no-merge-siblings has no effect in dense mode: sibling merging requires the sparse sequence representation"
      );
    }
    if !params.topology_ops.flip_parent_child {
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
    normalize_partition_rates(&mut [(length, reconstruction.gtr_mut())], &mut branch_lengths);
  }

  {
    let total_length = reconstruction.sequence_length();
    let indel_counts = gather_edge_indel_counts(&input.graph, &reconstruction);
    let sub_counts = gather_edge_sub_counts(&input.graph, &reconstruction).map_err(OperationError::classify)?;
    let effective_lengths =
      gather_edge_effective_lengths(&input.graph, &reconstruction).map_err(OperationError::classify)?;
    apply_initial_guess_mode(
      &input.graph,
      total_length,
      &indel_counts,
      &sub_counts,
      &effective_lengths,
      params.initial_guess,
      params.no_indels,
      &mut branch_lengths,
      names,
      log,
    )
    .map_err(OperationError::classify)?;
  }

  if let Some(spec) = &params.reroot_spec {
    progress_info!(log, "Rerooting before optimization: {spec:?}");
    stages.report("Rerooting", 0.2, "");
    reconstruction = pre_reroot_optimize(
      &input.graph,
      reconstruction,
      params.opt_method,
      params.no_indels,
      &mut branch_lengths,
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
    params.max_iter,
    params.dp,
    params.damping,
    params.opt_method,
    params.no_indels,
    params.topology_ops,
    branch_lengths,
    &loop_names,
  )
  .map_err(OperationError::classify)?;
  let branch_lengths = loop_result.branch_lengths;

  progress_info!(log, "Re-running marginal to populate subs_ml after optimization loop");
  let (reconstruction, _) = loop_result
    .reconstruction
    .marginal_update(&input.graph, &branch_lengths_or_zero(&branch_lengths))
    .map_err(OperationError::classify)?;

  Ok(OptimizeOutput {
    graph: input.graph,
    gtr: reconstruction.gtr().clone(),
    model_name: params.model,
    reconstruction,
    branch_lengths,
    names: loop_result.names,
  })
}

pub struct OptimizeParams {
  pub model: GtrModelName,
  pub dense: Option<bool>,
  pub max_iter: usize,
  pub dp: f64,
  pub damping: f64,
  pub opt_method: BranchOptMethod,
  pub initial_guess: InitialGuessMode,
  pub no_indels: bool,
  pub reroot_spec: Option<RerootSpec>,
  pub topology_ops: TopologyOps,
}

pub struct OptimizeInput {
  pub graph: Graph,
  pub alphabet: Alphabet,
  pub sequences: Vec<AlignmentRecord>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

#[derive(Debug, Serialize)]
pub struct OptimizeOutput {
  #[serde(skip)]
  pub graph: Graph,
  #[serde(skip)]
  pub gtr: GTR,
  pub model_name: GtrModelName,
  #[serde(skip)]
  pub reconstruction: MarginalReconstruction,
  #[serde(skip)]
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  #[serde(skip)]
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
}

fn pre_reroot_optimize(
  graph: &Graph,
  reconstruction: MarginalReconstruction,
  opt_method: BranchOptMethod,
  no_indels: bool,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<MarginalReconstruction, Report> {
  let old_branch_lengths = branch_lengths.clone();

  {
    let total_length = reconstruction.sequence_length();
    let indel_counts = gather_edge_indel_counts(graph, &reconstruction);
    if no_indels {
      run_optimize_mixed_inner(
        graph,
        total_length,
        &reconstruction,
        &indel_counts,
        opt_method,
        0.0,
        true,
        branch_lengths,
      )?;
    } else {
      run_optimize_mixed(
        graph,
        total_length,
        &reconstruction,
        &indel_counts,
        opt_method,
        branch_lengths,
      )?;
    }
  }

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
      let tip_keys = resolve_tip_keys(graph, tips, names)?;
      let mrca = common_ancestor(graph, &tip_keys)?;
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

fn resolve_tip_keys(
  graph: &Graph,
  tips: &[String],
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<Vec<GraphNodeKey>, Report> {
  tips
    .iter()
    .map(|tip| {
      names
        .iter()
        .find(|(_, name)| name.as_deref() == Some(tip.as_str()))
        .map(|(key, _)| *key)
        .ok_or_else(|| input_error(format!("Reroot tip not found: {tip}")))
    })
    .collect()
}
