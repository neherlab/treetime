use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::commands::optimize::args::{InitialGuessMode, TreetimeOptimizeArgs};
use crate::commands::optimize::optimize_unified::{initial_guess_mixed, run_optimize_mixed};
use crate::commands::optimize::partition_ops::{PartitionOptimizeOps, PartitionOptimizeVec};
use crate::commands::prune::run::merge_shared_mutation_branches;
use crate::gtr::get_gtr::{JC69Params, get_gtr_dense, get_gtr_sparse, jc69, write_gtr_json};
use crate::representation::algo::infer_dense::infer_dense;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition::traits::PartitionBranchOps;
use crate::representation::payload::ancestral::{GraphAncestral, annotate_branch_mutations};
use crate::seq::mutation::compose_substitutions;
use eyre::Report;
use itertools::{Itertools, chain, izip};
use log::debug;
use maplit::btreemap;
use num_traits::pow::pow;
use parking_lot::RwLock;
use serde::Serialize;
use std::path::Path;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_io::fasta::read_many_fasta;
use treetime_io::nex::{NexWriteOptions, nex_write_file};
use treetime_io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_read_file, nwk_write_file};
use treetime_utils::fmt::float::float_to_significant_digits;
use treetime_utils::make_error;

#[derive(Clone, Debug, Default)]
pub struct TreetimeOptimizeParams {
  pub sample_from_profile: bool,
  pub fixed_pi: bool,
}

pub fn run_optimize(args: &TreetimeOptimizeArgs) -> Result<(), Report> {
  let TreetimeOptimizeArgs {
    input_fastas,
    tree,
    alphabet,
    model_name,
    dense,
    outdir,
    max_iter,
    dp,
    damping,
    initial_guess,
  } = args;

  if !(0.0..1.0).contains(damping) {
    return make_error!("--damping must be in [0.0, 1.0), got {damping}");
  }

  let dense = dense.unwrap_or_else(infer_dense);
  let alphabet = Alphabet::new(alphabet.unwrap_or_default())?;
  let aln = read_many_fasta(input_fastas, &alphabet)?;
  let mut graph: GraphAncestral = nwk_read_file(tree)?;

  let sparse_partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>> = if !dense {
    #[allow(clippy::iter_on_single_items)]
    let partitions = [PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?, // FIXME: dummy temporary gtr should not be needed here
      alphabet: alphabet.clone(),
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }]
    .into_iter()
    .map(|p| Arc::new(RwLock::new(p)))
    .collect_vec();

    compress_sequences(&graph, &partitions, &aln)?;

    // FIXME: chicken & egg problem: to get a gtr we need partitions, to get partitions we need a gtr
    // FIXME: spaghetti code: dummy gtr is replaced by real gtr here
    for partition in &partitions {
      let gtr = get_gtr_sparse(model_name, partition, &graph)?;
      write_gtr_json(&gtr, *model_name, outdir, None)?;
      partition.write_arc().gtr = gtr;
    }

    partitions
  } else {
    vec![]
  };

  let dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>> = if dense {
    #[allow(clippy::iter_on_single_items)]
    let partitions = [PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?, // FIXME: dummy temporary gtr should not be needed here
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }]
    .into_iter()
    .map(|p| Arc::new(RwLock::new(p)))
    .collect_vec();

    partitions
  } else {
    vec![]
  };

  // Initialize marginal data for whichever partitions are active
  update_marginal(&graph, &sparse_partitions)?;
  if !dense_partitions.is_empty() {
    initialize_marginal(&graph, &dense_partitions, &aln)?;

    // FIXME: chicken & egg problem: to get a gtr we need partitions, to get partitions we need a gtr
    // Dense GTR requires profiles populated via initialize_marginal + update_marginal.
    // Run first marginal pass with dummy GTR, then replace with real model.
    update_marginal(&graph, &dense_partitions)?;
    for partition in &dense_partitions {
      let gtr = get_gtr_dense(model_name, partition, &graph)?;
      write_gtr_json(&gtr, *model_name, outdir, None)?;
      partition.write_arc().gtr = gtr;
    }

    // Re-run with inferred GTR so node posteriors reflect the real model.
    update_marginal(&graph, &dense_partitions)?;
  }

  let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);

  if should_run_initial_guess(*initial_guess, &graph) {
    initial_guess_mixed(&graph, &mixed_partitions)?;
  }

  let mut lh_prev = f64::MIN;
  for i in 0..*max_iter {
    let sparse_lh = update_marginal(&graph, &sparse_partitions)?;
    let dense_lh = update_marginal(&graph, &dense_partitions)?;
    let total_lh = sparse_lh + dense_lh;

    debug!(
      "Iteration {}: likelihood {} (sparse: {}, dense: {})",
      i + 1,
      float_to_significant_digits(total_lh, 7),
      float_to_significant_digits(sparse_lh, 7),
      float_to_significant_digits(dense_lh, 7)
    );

    if (total_lh - lh_prev).abs() < dp.abs() {
      break;
    }

    let old_branch_lengths = save_branch_lengths(&graph);
    run_optimize_mixed(&graph, &mixed_partitions)?;

    // Identify zero-optimal internal edges BEFORE damping blends them with old values
    let zero_optimal_edges = find_zero_optimal_internal_edges(&graph);

    apply_damping(&graph, &old_branch_lengths, *damping, i);

    // Collapse zero-optimal edges and merge shared mutations in resulting polytomies
    prune_and_merge_in_loop(&mut graph, &sparse_partitions, &dense_partitions, &zero_optimal_edges)?;

    lh_prev = total_lh;
  }

  // Annotate mutations from whichever partition type is active
  let branch_ops: Vec<Arc<RwLock<dyn PartitionBranchOps>>> = if dense {
    dense_partitions
      .iter()
      .cloned()
      .map(|p| -> Arc<RwLock<dyn PartitionBranchOps>> { p })
      .collect_vec()
  } else {
    sparse_partitions
      .iter()
      .cloned()
      .map(|p| -> Arc<RwLock<dyn PartitionBranchOps>> { p })
      .collect_vec()
  };
  annotate_branch_mutations(&graph, &branch_ops)?;

  write_graph(outdir, &graph)?;
  Ok(())
}

pub(crate) fn collect_optimize_partitions(
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
) -> PartitionOptimizeVec {
  chain!(
    dense_partitions
      .iter()
      .cloned()
      .map(|partition| -> Arc<RwLock<dyn PartitionOptimizeOps>> { partition }),
    sparse_partitions
      .iter()
      .cloned()
      .map(|partition| -> Arc<RwLock<dyn PartitionOptimizeOps>> { partition })
  )
  .collect_vec()
}

/// Save current branch lengths for all edges in graph traversal order.
pub(crate) fn save_branch_lengths<N, E, D>(graph: &Graph<N, E, D>) -> Vec<f64>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  graph
    .get_edges()
    .iter()
    .map(|edge_ref| edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
    .collect_vec()
}

/// Blend optimized branch lengths with saved old values using exponential damping.
///
/// At iteration `i` (0-based), each branch length becomes:
///   bl = bl_optimized * (1 - damping^(i+1)) + bl_old * damping^(i+1)
///
/// When `damping == 0.0`, damping_factor = 0 and the optimized value is kept unchanged.
/// Early iterations take conservative steps; later iterations approach the full update.
pub(crate) fn apply_damping<N, E, D>(graph: &Graph<N, E, D>, old_branch_lengths: &[f64], damping: f64, iteration: usize)
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  if damping == 0.0 {
    return;
  }
  debug_assert_eq!(
    graph.get_edges().len(),
    old_branch_lengths.len(),
    "apply_damping: edge count changed between save and apply"
  );
  let damping_factor = pow(damping, iteration + 1);
  let new_weight = 1.0 - damping_factor;
  for (edge_ref, &old_bl) in izip!(graph.get_edges(), old_branch_lengths) {
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let optimized_bl = edge.branch_length().unwrap_or(0.0);
    let damped_bl = optimized_bl * new_weight + old_bl * damping_factor;
    edge.set_branch_length(Some(damped_bl));
  }
}

/// Collect internal edge keys whose branch length is exactly zero after optimization.
///
/// These edges were identified as zero-optimal by `is_zero_branch_optimal()` during
/// `run_optimize_mixed()`. Call this BEFORE `apply_damping()`, which would blend the
/// zero values with old branch lengths and obscure the optimizer's decision.
///
/// Only internal (non-leaf-targeting) edges are candidates: collapsing a leaf edge
/// would remove the leaf from the tree.
pub(crate) fn find_zero_optimal_internal_edges(graph: &GraphAncestral) -> Vec<GraphEdgeKey> {
  graph
    .get_edges()
    .iter()
    .filter_map(|edge_ref| {
      let edge = edge_ref.read_arc();
      let bl = edge.payload().read_arc().branch_length().unwrap_or(f64::NAN);
      let target_is_leaf = graph.is_leaf(edge.target());
      (bl == 0.0 && !target_is_leaf).then_some(edge.key())
    })
    .collect_vec()
}

/// Collapse a single zero-optimal internal edge, updating both sparse and dense partition data.
///
/// Graph topology: the target node is removed and its children become children of the
/// source node. Edge keys for the former children are preserved (graph reuses them).
///
/// Sparse partitions: substitutions on the removed edge are composed with each child
/// edge's substitutions. Indels are concatenated (parent-first).
///
/// Dense partitions: stale node/edge entries are removed. Messages will be recomputed
/// by `update_marginal()` in the next iteration.
pub(crate) fn collapse_edge_for_optimize(
  graph: &mut GraphAncestral,
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  edge_key: GraphEdgeKey,
) -> Result<(), Report> {
  let target_node_key = graph.get_target_node_key(edge_key)?;

  let (_, removed_edge, new_edges) = graph.collapse_edge(edge_key)?;
  let removed_edge = removed_edge.payload().read_arc();

  for new_edge in &new_edges {
    let new_edge_key = new_edge.read_arc().key();
    let mut new_edge_payload = new_edge.write_arc().payload().write_arc();

    // Sum branch lengths: collapsed edge (0.0) + child edge = child edge unchanged
    if let (Some(bl1), Some(bl2)) = (removed_edge.branch_length(), new_edge_payload.branch_length()) {
      new_edge_payload.set_branch_length(Some(bl1 + bl2));
    }

    // Compose substitutions for sparse partitions
    for partition in sparse_partitions {
      let mut partition = partition.write_arc();
      let removed_data = partition.edges.get(&edge_key).cloned().unwrap_or_default();
      let child_edge = partition.edges.entry(new_edge_key).or_default();
      child_edge.subs = compose_substitutions(&removed_data.subs, &child_edge.subs)?;
      let mut merged_indels = removed_data.indels;
      merged_indels.append(&mut child_edge.indels);
      child_edge.indels = merged_indels;
    }
  }

  // Clean up stale partition data for the removed node and edge
  for partition in sparse_partitions {
    let mut partition = partition.write_arc();
    partition.nodes.remove(&target_node_key);
    partition.edges.remove(&edge_key);
  }
  for partition in dense_partitions {
    let mut partition = partition.write_arc();
    partition.nodes.remove(&target_node_key);
    partition.edges.remove(&edge_key);
  }

  Ok(())
}

/// Collapse zero-optimal internal edges and merge shared mutations in resulting polytomies.
///
/// Analogous to v0's `prune_short_branches()` inside the optimization loop:
/// edges whose optimal branch length is zero are collapsed, simplifying the tree
/// progressively across iterations. After collapsing, sibling branches in newly
/// formed polytomies that share identical substitutions are merged under new
/// internal nodes (sparse partitions only, requires discrete substitution data).
///
/// Returns true if any topology change occurred.
pub(crate) fn prune_and_merge_in_loop(
  graph: &mut GraphAncestral,
  sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  zero_optimal_edges: &[GraphEdgeKey],
) -> Result<bool, Report> {
  if zero_optimal_edges.is_empty() {
    return Ok(false);
  }

  // Override damped branch lengths back to zero for edges the optimizer identified
  // as zero-optimal. Damping is a convergence aid for continuous values; it should
  // not prevent collapsing degenerate edges.
  for &edge_key in zero_optimal_edges {
    if let Some(edge) = graph.get_edge(edge_key) {
      edge.write_arc().payload().write_arc().set_branch_length(Some(0.0));
    }
  }

  let mut collapsed = 0_usize;
  for &edge_key in zero_optimal_edges {
    // Edge may already be gone if a prior collapse in this batch removed it
    // (e.g., collapsing a parent also removed a child edge that was zero-optimal)
    if graph.get_edge(edge_key).is_none() {
      continue;
    }
    collapse_edge_for_optimize(graph, sparse_partitions, dense_partitions, edge_key)?;
    collapsed += 1;
  }

  if collapsed == 0 {
    return Ok(false);
  }

  debug!("Collapsed {collapsed} zero-optimal internal edges");

  // Merge shared mutations in newly formed polytomies (sparse only)
  if !sparse_partitions.is_empty() {
    merge_shared_mutation_branches(graph, sparse_partitions)?;
  }

  graph.build()?;
  Ok(true)
}

/// Decide whether to run the discrete-count initial guess based on mode and tree state.
fn should_run_initial_guess<N, E, D>(mode: InitialGuessMode, graph: &Graph<N, E, D>) -> bool
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  match mode {
    InitialGuessMode::Always => true,
    InitialGuessMode::Never => false,
    InitialGuessMode::Auto => graph
      .get_edges()
      .iter()
      .any(|edge_ref| edge_ref.read_arc().payload().read_arc().branch_length().is_none()),
  }
}

fn write_graph<N, E, D>(outdir: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk + Serialize,
  E: GraphEdge + EdgeToNwk + Serialize,
  D: Send + Sync + Default + Serialize,
{
  // json_write_file(
  //   outdir.as_ref().join("annotated_tree.graph.json"),
  //   &graph,
  //   JsonPretty(true),
  // )?;

  nwk_write_file(
    outdir.as_ref().join("annotated_tree.nwk"),
    graph,
    &NwkWriteOptions::default(),
  )?;

  nex_write_file(
    outdir.as_ref().join("annotated_tree.nexus"),
    graph,
    &NexWriteOptions::default(),
  )?;

  Ok(())
}
