use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::marginal::profile_branch_lengths;
use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
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
use treetime_io::fasta::FastaRecord;
use treetime_utils::{make_error, make_report};

/// Damping applied to the single pre-reroot branch-length pass, mirroring the
/// timetree pre-step: new lengths are blended `0.25 * new + 0.75 * old` so the
/// root search runs on stable distances without overshooting.
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
  /// Reroot policy before optimization. `None` keeps the input root.
  pub reroot_spec: Option<RerootSpec>,
  /// Which per-iteration topology-cleanup steps run. All on by default.
  pub topology_ops: TopologyOps,
}

pub struct OptimizeInput {
  pub graph: Graph,
  pub alphabet: Alphabet,
  pub sequences: Vec<FastaRecord>,
  /// Raw per-edge branch lengths captured from the Newick parse, keyed by edge id. The optimize loop
  /// takes ownership and makes it the source of truth: the initial guess, reroot, and per-edge
  /// optimizer all update it, and it exits as `OptimizeOutput.branch_lengths`.
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
  /// Final optimized branch lengths, keyed by edge id. The optimize loop is the source of truth;
  /// the command gather reads these.
  #[serde(skip)]
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  /// Final node names, keyed by node id, captured after the optimize loop's topology cleanup
  /// re-runs `assign_node_names`. The command gather and output writers read these.
  #[serde(skip)]
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
}

pub fn run(
  params: &OptimizeParams,
  mut input: OptimizeInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  progress: &dyn ProgressSink,
) -> Result<OptimizeOutput, Report> {
  if !(0.0..1.0).contains(&params.damping) {
    return make_error!("damping must be in [0.0, 1.0), got {}", params.damping);
  }

  let mut branch_lengths = std::mem::take(&mut input.branch_lengths);

  let created = create_marginal_partition(
    &input.graph,
    0,
    input.alphabet,
    &input.sequences,
    params.model,
    params.dense,
    &branch_lengths,
    names,
  )?;
  let model_name = created.model_name;

  let mut sparse_partitions: Vec<SparseReconstruction>;
  let mut dense_partitions: Vec<DenseReconstruction>;

  match created.partition {
    MarginalPartition::Sparse(partition, node_states) => {
      sparse_partitions = vec![SparseReconstruction::seeded(partition, node_states)];
      dense_partitions = vec![];
    },
    MarginalPartition::Dense(partition) => {
      // Dense leaf states are attached up front; the marginal passes then populate internal states.
      let node_states = partition.attach_sequences(&input.graph, &input.sequences, names)?;
      dense_partitions = vec![DenseReconstruction::seeded(partition, node_states)];
      sparse_partitions = vec![];
    },
  }

  // `merge_siblings` and `flip_parent_child` operate on discrete per-edge mutation lists, which
  // only sparse partitions carry. Under a dense build these steps never run, so disabling them
  // has no effect; warn the user their flag did nothing rather than failing silently.
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

  let profile_lengths = profile_branch_lengths(&branch_lengths);
  (sparse_partitions, _) = marginal_update_sparse(&input.graph, &profile_lengths, sparse_partitions)?;
  (dense_partitions, _) = marginal_update_dense(&input.graph, &profile_lengths, dense_partitions)?;

  if model_name == GtrModelName::Infer {
    let mut sparse_models: Vec<(usize, &mut GTR)> = sparse_partitions
      .iter_mut()
      .map(|family| (family.partition.length, family.partition.gtr_mut()))
      .collect();
    normalize_partition_rates(&mut sparse_models, &mut branch_lengths);
    let mut dense_models: Vec<(usize, &mut GTR)> = dense_partitions
      .iter_mut()
      .map(|family| (family.partition.length, family.partition.gtr_mut()))
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
    // Two-phase pattern (cf. timetree): optimize branch lengths first so the root
    // search uses stable distances, reroot, then re-optimize in the main loop
    // because the root move reshapes the optimization landscape.
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

  // Post-reroot names for the optimization loop. Reproject the threaded names onto the current node
  // set: a rerooted tree drops the old trivial root and adds an `N::default` split node (absent from
  // the map, so `None`); without a reroot the projection is the threaded map unchanged. The loop
  // refreshes this from each `assign_node_names` its topology cleanup runs.
  let loop_names: BTreeMap<GraphNodeKey, Option<String>> = input
    .graph
    .get_nodes()
    .iter()
    .map(|node| {
      let key = node.read_arc().key();
      (key, names.get(&key).cloned().flatten())
    })
    .collect();

  progress.check_cancelled()?;
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
  let marginal_bl = profile_branch_lengths(&branch_lengths);
  let (sparse_partitions, _) = marginal_update_sparse(&input.graph, &marginal_bl, loop_result.sparse_partitions)?;
  let (dense_partitions, _) = marginal_update_dense(&input.graph, &marginal_bl, loop_result.dense_partitions)?;

  // Read the GTR back from the owning partition, not from a snapshot taken at
  // creation time. For `--gtr=infer` the partition's `mu` is normalized to 1.0
  // by `normalize_partition_rates` above (rate absorbed into branch lengths);
  // a creation-time clone would still carry the raw inferred `mu` and disagree
  // with the rate-scaled branch lengths shipped alongside it.
  let gtr = if let Some(family) = sparse_partitions.first() {
    family.partition.gtr().clone()
  } else if let Some(family) = dense_partitions.first() {
    family.partition.gtr().clone()
  } else {
    return make_error!("optimize produced no partition to read the GTR from");
  };

  // The loop's topology cleanup collapses edges, resolves polytomies, and re-runs `assign_node_names`
  // (adding or removing node keys and naming new internal nodes); the loop refreshes and returns the
  // node-name map from each such call, so the command gather and output writers read the final tree's
  // names.
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

/// Single damped branch-length pass run before rerooting.
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
  let profile_lengths = profile_branch_lengths(branch_lengths);
  let (sparse_partitions, _) = marginal_update_sparse(graph, &profile_lengths, sparse_partitions)?;
  let (dense_partitions, _) = marginal_update_dense(graph, &profile_lengths, dense_partitions)?;
  Ok((sparse_partitions, dense_partitions))
}

/// Reroot the tree by the requested date-free policy and reconcile partitions.
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
    .map(|family| reroot_sparse(family.partition, family.node_states, &changes))
    .try_collect()?;
  let dense_partitions: Vec<_> = dense_partitions
    .into_iter()
    .map(|family| reroot_dense(family.partition, family.node_states, &changes))
    .collect();

  let profile_lengths = profile_branch_lengths(branch_lengths);
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
