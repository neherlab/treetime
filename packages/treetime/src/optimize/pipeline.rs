use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::optimize::dispatch::{run_optimize_mixed, run_optimize_mixed_inner};
use crate::optimize::iteration::{apply_damping, save_branch_lengths};
use crate::optimize::params::{BranchOptMethod, InitialGuessMode};
use crate::optimize::run_loop::{
  apply_initial_guess_mode, collect_optimize_partitions, normalize_partition_rates, run_optimize_loop,
};
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::partition::traits::{HasGtr, PartitionOptimizeVec, PartitionRerootOps};
use crate::payload::ancestral::GraphAncestral;
use crate::progress::ProgressSink;
use crate::reroot::div_stats::DivStats;
use crate::reroot::div_stats_traversal::compute_div_stats;
use crate::reroot::orchestrate::{RerootTopologyParams, reroot_at_node, reroot_in_place};
use crate::reroot::params::BrentParams;
use crate::reroot::variance::VarianceModel;
use eyre::Report;
use log::info;
use parking_lot::RwLock;
use serde::Serialize;
use std::sync::Arc;
use treetime_graph::common_ancestor::common_ancestor;
use treetime_graph::node::{GraphNodeKey, Named};
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
}

pub struct OptimizeInput {
  pub graph: GraphAncestral,
  pub alphabet: Alphabet,
  pub sequences: Vec<FastaRecord>,
}

#[derive(Debug, Serialize)]
pub struct OptimizeOutput {
  #[serde(skip)]
  pub graph: GraphAncestral,
  #[serde(skip)]
  pub gtr: GTR,
  pub model_name: GtrModelName,
  #[serde(skip)]
  pub sparse_partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>>,
  #[serde(skip)]
  pub dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>>,
}

pub fn run(
  params: &OptimizeParams,
  mut input: OptimizeInput,
  progress: &dyn ProgressSink,
) -> Result<OptimizeOutput, Report> {
  if !(0.0..1.0).contains(&params.damping) {
    return make_error!("damping must be in [0.0, 1.0), got {}", params.damping);
  }

  let created = create_marginal_partition(
    &input.graph,
    0,
    input.alphabet,
    &input.sequences,
    params.model,
    params.dense,
  )?;
  let model_name = created.model_name;

  let sparse_partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>>;
  let dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>>;

  match created.partition {
    MarginalPartition::Sparse(p) => {
      sparse_partitions = vec![Arc::new(RwLock::new(p))];
      dense_partitions = vec![];
    },
    MarginalPartition::Dense(p) => {
      sparse_partitions = vec![];
      dense_partitions = vec![Arc::new(RwLock::new(p))];
    },
  }

  update_marginal(&input.graph, &sparse_partitions)?;
  if !dense_partitions.is_empty() {
    initialize_marginal(&input.graph, &dense_partitions, &input.sequences)?;
    update_marginal(&input.graph, &dense_partitions)?;
  }

  let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);

  if model_name == GtrModelName::Infer {
    normalize_partition_rates(&input.graph, &sparse_partitions);
    normalize_partition_rates(&input.graph, &dense_partitions);
  }

  apply_initial_guess_mode(&input.graph, &mixed_partitions, params.initial_guess, params.no_indels)?;

  if let Some(spec) = &params.reroot_spec {
    info!("Rerooting before optimization: {spec:?}");
    progress.report("Rerooting", 0.2, "");
    // Two-phase pattern (cf. timetree): optimize branch lengths first so the root
    // search uses stable distances, reroot, then re-optimize in the main loop
    // because the root move reshapes the optimization landscape.
    pre_reroot_optimize(
      &input.graph,
      &mixed_partitions,
      &sparse_partitions,
      &dense_partitions,
      params.opt_method,
      params.no_indels,
    )?;
    reroot_optimize(&mut input.graph, spec, &sparse_partitions, &dense_partitions)?;
  }

  progress.check_cancelled()?;
  progress.report("Optimizing branch lengths", 0.3, "");
  run_optimize_loop(
    &mut input.graph,
    &sparse_partitions,
    &dense_partitions,
    &mixed_partitions,
    params.max_iter,
    params.dp,
    params.damping,
    params.opt_method,
    params.no_indels,
  )?;

  info!("Re-running marginal to populate subs_ml after optimization loop");
  update_marginal(&input.graph, &sparse_partitions)?;
  update_marginal(&input.graph, &dense_partitions)?;

  // Read the GTR back from the owning partition, not from a snapshot taken at
  // creation time. For `--gtr=infer` the partition's `mu` is normalized to 1.0
  // by `normalize_partition_rates` above (rate absorbed into branch lengths);
  // a creation-time clone would still carry the raw inferred `mu` and disagree
  // with the rate-scaled branch lengths shipped alongside it.
  let gtr = if let Some(p) = sparse_partitions.first() {
    p.read_arc().gtr().clone()
  } else if let Some(p) = dense_partitions.first() {
    p.read_arc().gtr().clone()
  } else {
    return make_error!("optimize produced no partition to read the GTR from");
  };

  Ok(OptimizeOutput {
    graph: input.graph,
    gtr,
    model_name,
    sparse_partitions,
    dense_partitions,
  })
}

/// Single damped branch-length pass run before rerooting.
fn pre_reroot_optimize(
  graph: &GraphAncestral,
  mixed_partitions: &PartitionOptimizeVec,
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  opt_method: BranchOptMethod,
  no_indels: bool,
) -> Result<(), Report> {
  let old_branch_lengths = save_branch_lengths(graph);

  if no_indels {
    run_optimize_mixed_inner(graph, mixed_partitions, opt_method, 0.0, true)?;
  } else {
    run_optimize_mixed(graph, mixed_partitions, opt_method)?;
  }

  apply_damping(graph, &old_branch_lengths, PRE_REROOT_DAMPING, 0);
  update_marginal(graph, sparse_partitions)?;
  update_marginal(graph, dense_partitions)?;
  Ok(())
}

/// Reroot the tree by the requested date-free policy and reconcile partitions.
fn reroot_optimize(
  graph: &mut GraphAncestral,
  spec: &RerootSpec,
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
) -> Result<(), Report> {
  let variance = VarianceModel::default();
  let topo = RerootTopologyParams::default();
  let opt_params = BrentParams::default();

  let reroot_result = match spec {
    RerootSpec::Method(RerootMethod::MinDev) => {
      let field = compute_div_stats(graph, &variance)?;
      reroot_in_place::<_, _, _, DivStats, _>(
        graph,
        &field.edge_stats,
        &field.root_stats,
        &variance,
        &opt_params,
        topo,
        |_graph, _inverted| Ok(()),
      )?
    },
    RerootSpec::Tips(tips) => {
      let tip_keys = resolve_tip_keys(graph, tips)?;
      let mrca = common_ancestor(graph, &tip_keys)?;
      reroot_at_node(graph, mrca, topo, |_graph, _inverted| Ok(()))?
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

  for partition in sparse_partitions {
    partition.write_arc().apply_reroot(&changes)?;
  }
  for partition in dense_partitions {
    partition.write_arc().apply_reroot(&changes)?;
  }

  update_marginal(graph, sparse_partitions)?;
  update_marginal(graph, dense_partitions)?;
  Ok(())
}

fn resolve_tip_keys(graph: &GraphAncestral, tips: &[String]) -> Result<Vec<GraphNodeKey>, Report> {
  tips
    .iter()
    .map(|tip| {
      graph
        .find_node(|node| node.name().is_some_and(|name| name.as_ref() == tip))
        .ok_or_else(|| make_report!("Reroot tip not found: {tip}"))
    })
    .collect()
}
