use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::optimize::params::{BranchOptMethod, InitialGuessMode};
use crate::optimize::run_loop::{
  apply_initial_guess_mode, collect_optimize_partitions, normalize_partition_rates, run_optimize_loop,
};
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::partition::traits::HasGtr;
use crate::payload::ancestral::GraphAncestral;
use crate::progress::ProgressSink;
use eyre::Report;
use log::info;
use parking_lot::RwLock;
use serde::Serialize;
use std::sync::Arc;
use treetime_io::fasta::FastaRecord;
use treetime_utils::make_error;

pub struct OptimizeParams {
  pub model: GtrModelName,
  pub dense: Option<bool>,
  pub max_iter: usize,
  pub dp: f64,
  pub damping: f64,
  pub opt_method: BranchOptMethod,
  pub initial_guess: InitialGuessMode,
  pub no_indels: bool,
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
