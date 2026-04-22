#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, load_alignment_for_dataset, load_dates_for_dataset,
  };
  use crate::commands::ancestral::fitch::compress_sequences;
  use crate::commands::ancestral::marginal::initialize_marginal;
  use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
  use crate::commands::clock::date_constraints::load_date_constraints;
  use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::commands::optimize::args::BranchOptMethod;
  use crate::commands::optimize::optimize_unified::run_optimize_mixed;
  use crate::commands::optimize::partition_ops::PartitionOptimizeOps;
  use crate::commands::timetree::inference::runner::run_timetree;
  use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
  use crate::commands::timetree::utils::{
    extract_node_times, initialize_clock_totals_from_time_distributions, initialize_node_divergences,
  };
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::representation::payload::timetree::{EdgeTimetree, NodeTimetree};
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use rstest::rstest;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::seq;

  fn extract_branch_lengths(graph: &GraphTimetree) -> Vec<f64> {
    graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .collect_vec()
  }

  /// Verify that the ML branch-length pre-step modifies branch lengths.
  ///
  /// Sets up a sparse partition with JC69, runs marginal reconstruction,
  /// then runs one pass of Brent optimization. The test asserts that at
  /// least some branch lengths differ from their input values.
  #[rustfmt::skip]
  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  #[trace]
  fn test_gm_runner_pre_optimize_changes_branch_lengths(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];

    let graph: GraphTimetree = nwk_read_str(case.rerooted_tree_nwk())?;
    let aln = load_alignment_for_dataset(dataset)?;
    let sparse_partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: ALPHABET.clone(),
      length: case.sequence_length(),
      nodes: btreemap! {},
      edges: btreemap! {},
      root_sequence: seq![],
    }));

    compress_sequences(&graph, std::slice::from_ref(&sparse_partition), &aln)?;
    sparse_partition.write_arc().extract_root_sequence(&graph)?;

    let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>> =
      vec![sparse_partition];
    initialize_marginal(&graph, &partitions, &aln)?;

    let before = extract_branch_lengths(&graph);

    // Run one pass of ML optimization (matching v0's optimize_tree(max_iter=1))
    #[allow(trivial_casts)]
    let opt_partitions: Vec<Arc<RwLock<dyn PartitionOptimizeOps>>> = partitions
      .iter()
      .map(|p| Arc::clone(p) as Arc<RwLock<dyn PartitionOptimizeOps>>)
      .collect();
    run_optimize_mixed(&graph, &opt_partitions, BranchOptMethod::BrentSqrt)?;

    let after = extract_branch_lengths(&graph);

    // At least some branch lengths should change
    let n_changed = before
      .iter()
      .zip(&after)
      .filter(|(b, a)| (*b - *a).abs() > 1e-10)
      .count();

    assert!(
      n_changed > 0,
      "Expected at least one branch length to change after ML optimization, but none did"
    );

    Ok(())
  }

  /// Verify that the full timetree pipeline succeeds with the pre-optimization step.
  ///
  /// Runs the same pipeline as marginal sparse but with the pre-optimization step
  /// inserted before time inference. Checks that the pipeline completes without
  /// error and produces node times.
  #[rustfmt::skip]
  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  #[trace]
  fn test_gm_runner_pre_optimize_pipeline_succeeds(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];

    let mut graph: GraphTimetree = nwk_read_str(case.rerooted_tree_nwk())?;
    let dates = load_dates_for_dataset(dataset)?;
    load_date_constraints(&dates, &graph)?;

    let aln = load_alignment_for_dataset(dataset)?;
    let sparse_partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: ALPHABET.clone(),
      length: case.sequence_length(),
      nodes: btreemap! {},
      edges: btreemap! {},
      root_sequence: seq![],
    }));

    compress_sequences(&graph, std::slice::from_ref(&sparse_partition), &aln)?;
    sparse_partition.write_arc().extract_root_sequence(&graph)?;

    let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>> =
      vec![sparse_partition];
    initialize_marginal(&graph, &partitions, &aln)?;
    initialize_node_divergences(&graph);

    // Pre-optimization step (matching v0 flow)
    #[allow(trivial_casts)]
    let opt_partitions: Vec<Arc<RwLock<dyn PartitionOptimizeOps>>> = partitions
      .iter()
      .map(|p| Arc::clone(p) as Arc<RwLock<dyn PartitionOptimizeOps>>)
      .collect();
    run_optimize_mixed(&graph, &opt_partitions, BranchOptMethod::BrentSqrt)?;
    crate::commands::ancestral::marginal::update_marginal(&graph, &partitions)?;

    let clock_model = estimate_clock_model_with_reroot(
      &mut graph,
      &ClockParams::default(),
      Some(case.clock_rate()),
      true,
      &BranchPointOptimizationParams::default(),
      None,
    )?;

    initialize_clock_totals_from_time_distributions(&graph)?;
    run_timetree(&mut graph, &partitions, &clock_model, None)?;

    let actual = extract_node_times(&graph);
    assert!(
      !actual.is_empty(),
      "Expected node times to be populated after timetree inference with pre-optimization"
    );

    // All internal nodes should have finite times
    for (name, time) in &actual {
      assert!(
        time.is_finite(),
        "Node {name} has non-finite time {time} after timetree inference with pre-optimization"
      );
    }

    Ok(())
  }
}
