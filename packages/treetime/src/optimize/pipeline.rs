use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::marginal::branch_lengths_or_zero;
use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
use crate::cancel::Cancel;
use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::optimize::dispatch::{run_optimize_mixed, run_optimize_mixed_inner};
use crate::optimize::gather::{
  gather_edge_contributions, gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts,
  total_sequence_length,
};
use crate::optimize::iteration::apply_damping;
use crate::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};
use crate::optimize::run_loop::{
  apply_initial_guess_mode, marginal_update_dense, marginal_update_sparse, normalize_partition_rates, run_optimize_loop,
};
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::marginal::dense::reroot::reroot_dense;
use crate::partition::marginal::sparse::reroot::reroot_sparse;
use crate::progress::ProgressSink;
use crate::reroot::div_stats::DivStats;
use crate::reroot::div_stats_traversal::compute_div_stats;
use crate::reroot::orchestrate::{RerootTopologyParams, reroot_at_node, reroot_in_place};
use crate::reroot::params::BrentParams;
use crate::reroot::variance::VarianceModel;
use crate::seq::alignment::node_seq_inputs;
use eyre::Report;
use itertools::Itertools;
use log::{info, warn};
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::common_ancestor::common_ancestor;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::RerootChanges;
use treetime_primitives::AlignmentRecord;
use treetime_utils::{make_error, make_report};

const PRE_REROOT_DAMPING: f64 = 0.75;

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
  pub sparse_partitions: Vec<SparseReconstruction>,
  #[serde(skip)]
  pub dense_partitions: Vec<DenseReconstruction>,
  #[serde(skip)]
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  #[serde(skip)]
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
}

pub fn run(
  params: &OptimizeParams,
  mut input: OptimizeInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<OptimizeOutput, OperationError> {
  if !(0.0..1.0).contains(&params.damping) {
    return Err(OperationError::InvalidParams(make_report!(
      "damping must be in [0.0, 1.0), got {}",
      params.damping
    )));
  }

  let mut branch_lengths = std::mem::take(&mut input.branch_lengths);
  let sequences = std::mem::take(&mut input.sequences);
  let node_inputs = node_seq_inputs(&input.graph, names, sequences);

  let created = create_marginal_partition(
    &input.graph,
    0,
    input.alphabet,
    &node_inputs,
    params.model,
    params.dense,
    &branch_lengths_or_zero(&branch_lengths),
  )?;
  let model_name = created.model_name;
  let gtr = created.gtr;

  let mut sparse_partitions: Vec<SparseReconstruction>;
  let mut dense_partitions: Vec<DenseReconstruction>;

  match created.partition {
    MarginalPartition::Sparse(partition, node_states) => {
      sparse_partitions = vec![SparseReconstruction::seeded(partition, gtr, node_states)];
      dense_partitions = vec![];
    },
    MarginalPartition::Dense(partition) => {
      let node_states = partition.attach_sequences(&input.graph, &node_inputs)?;
      dense_partitions = vec![DenseReconstruction::seeded(partition, gtr, node_states)];
      sparse_partitions = vec![];
    },
  }

  if sparse_partitions.is_empty() {
    if !params.topology_ops.merge_siblings {
      warn!(
        "--no-merge-siblings has no effect in dense mode: sibling merging requires the sparse sequence representation"
      );
    }
    if !params.topology_ops.flip_parent_child {
      warn!(
        "--no-flip-parent-child has no effect in dense mode: the reversion hoist requires the sparse sequence representation"
      );
    }
  }

  let profile_lengths = branch_lengths_or_zero(&branch_lengths);
  (sparse_partitions, _) = marginal_update_sparse(&input.graph, &profile_lengths, sparse_partitions)?;
  (dense_partitions, _) = marginal_update_dense(&input.graph, &profile_lengths, dense_partitions)?;

  if model_name == GtrModelName::Infer {
    let mut sparse_models: Vec<(usize, &mut GTR)> = sparse_partitions
      .iter_mut()
      .map(|family| (family.partition.length, &mut family.gtr))
      .collect();
    normalize_partition_rates(&mut sparse_models, &mut branch_lengths);
    let mut dense_models: Vec<(usize, &mut GTR)> = dense_partitions
      .iter_mut()
      .map(|family| (family.partition.length, &mut family.gtr))
      .collect();
    normalize_partition_rates(&mut dense_models, &mut branch_lengths);
  }

  {
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let indel_counts = gather_edge_indel_counts(&input.graph, &dense_partitions, &sparse_partitions);
    let sub_counts = gather_edge_sub_counts(&input.graph, &dense_partitions, &sparse_partitions)?;
    let effective_lengths = gather_edge_effective_lengths(&input.graph, &dense_partitions, &sparse_partitions)?;
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
    )?;
  }

  if let Some(spec) = &params.reroot_spec {
    info!("Rerooting before optimization: {spec:?}");
    progress.report("Rerooting", 0.2, "");
    (sparse_partitions, dense_partitions) = pre_reroot_optimize(
      &input.graph,
      sparse_partitions,
      dense_partitions,
      params.opt_method,
      params.no_indels,
      &mut branch_lengths,
    )?;
    (sparse_partitions, dense_partitions) = reroot_optimize(
      &mut input.graph,
      spec,
      sparse_partitions,
      dense_partitions,
      &mut branch_lengths,
      names,
    )?;
  }

  let loop_names: BTreeMap<GraphNodeKey, Option<String>> = input
    .graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      (key, names.get(&key).cloned().flatten())
    })
    .collect();

  cancel.check()?;
  progress.report("Optimizing branch lengths", 0.3, "");
  let loop_result = run_optimize_loop(
    &mut input.graph,
    sparse_partitions,
    dense_partitions,
    params.max_iter,
    params.dp,
    params.damping,
    params.opt_method,
    params.no_indels,
    params.topology_ops,
    branch_lengths,
    &loop_names,
  )?;
  let branch_lengths = loop_result.branch_lengths;

  info!("Re-running marginal to populate subs_ml after optimization loop");
  let marginal_bl = branch_lengths_or_zero(&branch_lengths);
  let (sparse_partitions, _) = marginal_update_sparse(&input.graph, &marginal_bl, loop_result.sparse_partitions)?;
  let (dense_partitions, _) = marginal_update_dense(&input.graph, &marginal_bl, loop_result.dense_partitions)?;

  let gtr = if let Some(family) = sparse_partitions.first() {
    family.gtr.clone()
  } else if let Some(family) = dense_partitions.first() {
    family.gtr.clone()
  } else {
    return Err(OperationError::InferenceFailed(make_report!(
      "optimize produced no partition to read the GTR from"
    )));
  };

  let names = loop_result.names;

  Ok(OptimizeOutput {
    graph: input.graph,
    gtr,
    model_name,
    sparse_partitions,
    dense_partitions,
    branch_lengths,
    names,
  })
}

fn pre_reroot_optimize(
  graph: &Graph,
  sparse_partitions: Vec<SparseReconstruction>,
  dense_partitions: Vec<DenseReconstruction>,
  opt_method: BranchOptMethod,
  no_indels: bool,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(Vec<SparseReconstruction>, Vec<DenseReconstruction>), Report> {
  let old_branch_lengths = branch_lengths.clone();

  {
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let contributions = gather_edge_contributions(graph, &dense_partitions, &sparse_partitions)?;
    let indel_counts = gather_edge_indel_counts(graph, &dense_partitions, &sparse_partitions);
    if no_indels {
      run_optimize_mixed_inner(
        graph,
        total_length,
        &contributions,
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
        &contributions,
        &indel_counts,
        opt_method,
        branch_lengths,
      )?;
    }
  }

  apply_damping(branch_lengths, &old_branch_lengths, PRE_REROOT_DAMPING, 0);
  let profile_lengths = branch_lengths_or_zero(branch_lengths);
  let (sparse_partitions, _) = marginal_update_sparse(graph, &profile_lengths, sparse_partitions)?;
  let (dense_partitions, _) = marginal_update_dense(graph, &profile_lengths, dense_partitions)?;
  Ok((sparse_partitions, dense_partitions))
}

fn reroot_optimize(
  graph: &mut Graph,
  spec: &RerootSpec,
  sparse_partitions: Vec<SparseReconstruction>,
  dense_partitions: Vec<DenseReconstruction>,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<(Vec<SparseReconstruction>, Vec<DenseReconstruction>), Report> {
  let variance = VarianceModel::default();
  let topo = RerootTopologyParams::default();
  let opt_params = BrentParams::default();

  let reroot_result = match spec {
    RerootSpec::Method(RerootMethod::MinDev) => {
      let field = compute_div_stats(graph, branch_lengths, &variance)?;
      reroot_in_place::<DivStats, _>(
        graph,
        &field.edge_stats,
        &field.root_stats,
        &variance,
        &opt_params,
        topo,
        branch_lengths,
        |_graph, _inverted| Ok(()),
      )?
    },
    RerootSpec::Tips(tips) => {
      let tip_keys = resolve_tip_keys(graph, tips, names)?;
      let mrca = common_ancestor(graph, &tip_keys)?;
      reroot_at_node(graph, mrca, topo, branch_lengths, |_graph, _inverted| Ok(()))?
    },
    RerootSpec::Method(method) => {
      return make_error!("optimize cannot reroot with a date-dependent method: {method:?}");
    },
  };

  let changes = RerootChanges {
    edge_split: reroot_result.edge_split,
    edge_merge: reroot_result.edge_merge,
    inverted_edge_keys: reroot_result.inverted_edge_keys,
  };

  let sparse_partitions: Vec<_> = sparse_partitions
    .into_iter()
    .map(|family| reroot_sparse(family.partition, family.gtr, family.node_states, &changes))
    .try_collect()?;
  let dense_partitions: Vec<_> = dense_partitions
    .into_iter()
    .map(|family| reroot_dense(family.partition, family.gtr, family.node_states, &changes))
    .collect();

  let profile_lengths = branch_lengths_or_zero(branch_lengths);
  let (sparse_partitions, _) = marginal_update_sparse(graph, &profile_lengths, sparse_partitions)?;
  let (dense_partitions, _) = marginal_update_dense(graph, &profile_lengths, dense_partitions)?;
  Ok((sparse_partitions, dense_partitions))
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
        .ok_or_else(|| make_report!("Reroot tip not found: {tip}"))
    })
    .collect()
}
