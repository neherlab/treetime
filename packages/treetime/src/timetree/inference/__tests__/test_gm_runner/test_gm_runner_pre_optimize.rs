#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, load_alignment_for_dataset, load_dates_for_dataset,
  };

  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::partition::timetree::marginal::{initialize_marginal_timetree, marginal_update_timetree};
  use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
  use crate::clock::clock_state::ClockState;
  use crate::clock::date_constraints::load_date_constraints;
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::clock::reroot::RerootParams;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::params::BranchOptMethod;
  use crate::partition::traits::PartitionOptimizeOps;
  use crate::timetree::inference::runner::run_timetree;
  use crate::timetree::timetree_state::TimetreeState;
  use crate::timetree::utils::{extract_node_times, initialize_node_divergences};
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;

  use crate::partition::timetree::partition::PartitionTimetree;
  use eyre::Report;
  use itertools::Itertools;
  use treetime_graph::graph::Graph;

  use rstest::rstest;

  use treetime_io::nwk::{NwkParse, nwk_read_str};

  fn extract_branch_lengths(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Vec<f64> {
    graph
      .get_edges()
      .iter()
      .map(|e| branch_lengths[&e.read_arc().key()].unwrap_or(0.0))
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

    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str(case.rerooted_tree_nwk())?;

    let graph: Graph = graph;
    let aln = load_alignment_for_dataset(dataset)?;
    let fitch = create_fitch_partition(&graph, 0, ALPHABET.clone(), &aln, &names)?;
    let (partition, node_states) = fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?;
    let sparse_partition = PartitionTimetree::Sparse(SparseReconstruction {
      partition,
      node_states,
      backward: BTreeMap::new(),
      forward: BTreeMap::new(),
      estimates: BTreeMap::new(),
    });

    let mut partitions: Vec<PartitionTimetree> = vec![sparse_partition];
    initialize_marginal_timetree(&graph, &profile_branch_lengths(&branch_lengths), &mut partitions, &aln, &names)?.value();

    let before = extract_branch_lengths(&graph, &branch_lengths);

    // Run one pass of ML optimization (matching v0's optimize_tree(max_iter=1))
    #[allow(trivial_casts)]
    let opt_partitions: Vec<&dyn PartitionOptimizeOps> = partitions.iter().map(|p| p as &dyn PartitionOptimizeOps).collect();
    run_optimize_mixed(&graph, &opt_partitions, BranchOptMethod::BrentSqrt, &mut branch_lengths)?;

    let after = extract_branch_lengths(&graph, &branch_lengths);

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

    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str(case.rerooted_tree_nwk())?;

    let mut graph: Graph = graph;
    let dates = load_dates_for_dataset(dataset)?;
    let constraints = load_date_constraints(&dates, &graph, &names)?;

    let aln = load_alignment_for_dataset(dataset)?;
    let fitch = create_fitch_partition(&graph, 0, ALPHABET.clone(), &aln, &names)?;
    let (partition, node_states) = fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?;
    let sparse_partition = PartitionTimetree::Sparse(SparseReconstruction {
      partition,
      node_states,
      backward: BTreeMap::new(),
      forward: BTreeMap::new(),
      estimates: BTreeMap::new(),
    });

    let mut partitions: Vec<PartitionTimetree> = vec![sparse_partition];
    initialize_marginal_timetree(&graph, &profile_branch_lengths(&branch_lengths), &mut partitions, &aln, &names)?.value();
    let mut clock_state = ClockState::new(&graph);
    initialize_node_divergences(&graph, &mut clock_state, &branch_lengths, &names)?;

    // Pre-optimization step (matching v0 flow)
    #[allow(trivial_casts)]
    let opt_partitions: Vec<&dyn PartitionOptimizeOps> = partitions.iter().map(|p| p as &dyn PartitionOptimizeOps).collect();
    run_optimize_mixed(&graph, &opt_partitions, BranchOptMethod::BrentSqrt, &mut branch_lengths)?;
    marginal_update_timetree(&graph, &profile_branch_lengths(&branch_lengths), &mut partitions)?;

    let times = TimetreeState::seed_from_values(&graph, &constraints).likely_times();
    let mut clock_estimate_state = ClockState::seed_from_values(&graph, &times);
    let names_tt_1 = names.clone();
    let clock_model = estimate_clock_model_with_reroot_policy(
      &mut graph,
      &mut clock_estimate_state,
      &ClockParams::default(),
      Some(case.clock_rate()),
      true,
      &BranchPointOptimizationParams::default(),
      &RerootParams::default(),
      &mut branch_lengths,
      None, &names_tt_1
    )?
    .into_clock_model()?;

    let mut state = TimetreeState::seed_from_values(&graph, &constraints);
    let run_branch_lengths = branch_lengths;
    let run_names = names.clone();
    run_timetree(
      &mut graph,
      &partitions,      &run_branch_lengths,
      &run_names,
      &clock_model,
      None,
      false,
      &mut state,
      &mut clock_state,
    )?;

    let actual = extract_node_times(&graph, &names, &state);
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
