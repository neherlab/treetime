#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
pub(super) mod tests {
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};

  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_partitions, simple_alignment,
  };
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::gather::{gather_edge_contributions, gather_edge_indel_counts, total_sequence_length};
  use crate::optimize::indel::{estimate_indel_rate, poisson_indel_log_lh};
  use crate::optimize::likelihood::{OptimizationMetrics, evaluate_mixed, evaluate_mixed_log_lh_only};

  use crate::optimize::method_newton::newton_tolerance_t;
  use crate::optimize::method_newton::{chain_rule_log, chain_rule_sqrt};
  use crate::optimize::params::BranchOptMethod;

  use crate::seq::indel::InDel;

  use eyre::Report;
  use helpers::*;

  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;

  use proptest::prelude::*;

  #[rustfmt::skip]
  #[rstest]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[case::brent(      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethod::BrentLog)]
  #[trace]
  fn test_optimize_method_equivalence_no_indels(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let aln = simple_alignment()?;
    let (dense_mixed_partitions, sparse_mixed_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;
    let total_length = total_sequence_length(&dense_mixed_partitions, &sparse_mixed_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_mixed_partitions, &sparse_mixed_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_mixed_partitions, &sparse_mixed_partitions);

    run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;

    for (i, edge_ref) in graph.get_edges().enumerate() {
      let bl = branch_lengths[&edge_ref.key()].unwrap_or(f64::NAN);
      assert!(bl.is_finite(), "Edge {i}: branch length is not finite ({bl})");
      assert!(bl >= 0.0, "Edge {i}: branch length is negative ({bl})");
    }

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[case::brent(      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethod::BrentLog)]
  #[trace]
  fn test_optimize_method_local_optimality(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (dense_mixed_partitions, sparse_mixed_partitions, indel_rate) = setup_with_indels(&graph, &names, &mut branch_lengths, 4)?;
    let total_length = total_sequence_length(&dense_mixed_partitions, &sparse_mixed_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_mixed_partitions, &sparse_mixed_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_mixed_partitions, &sparse_mixed_partitions);

    run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;

    let bl = first_edge_bl(&graph, &branch_lengths);
    assert!(bl > 0.0 && bl.is_finite(), "Optimized BL must be positive and finite, got {bl}");

    let lh_opt = eval_combined_first_edge(&graph, &dense_mixed_partitions, &sparse_mixed_partitions, indel_rate, bl)?;

    for &frac in &[0.001, 0.01, 0.1] {
      let delta = bl * frac;
      if bl - delta > 0.0 {
        let lh_below = eval_combined_first_edge(&graph, &dense_mixed_partitions, &sparse_mixed_partitions, indel_rate, bl - delta)?;
        assert!(
          lh_opt >= lh_below - 1e-10,
          "{method:?}: lh at t*={bl} ({lh_opt}) < lh at t*-{frac}*t ({lh_below})"
        );
      }
      let lh_above = eval_combined_first_edge(&graph, &dense_mixed_partitions, &sparse_mixed_partitions, indel_rate, bl + delta)?;
      assert!(
        lh_opt >= lh_above - 1e-10,
        "{method:?}: lh at t*={bl} ({lh_opt}) < lh at t*+{frac}*t ({lh_above})"
      );
    }

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[trace]
  fn test_optimize_method_stationarity(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (dense_mixed_partitions, sparse_mixed_partitions, indel_rate) = setup_with_indels(&graph, &names, &mut branch_lengths, 2)?;
    let total_length = total_sequence_length(&dense_mixed_partitions, &sparse_mixed_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_mixed_partitions, &sparse_mixed_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_mixed_partitions, &sparse_mixed_partitions);

    run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;

    let bl = first_edge_bl(&graph, &branch_lengths);
    let metrics = eval_metrics_first_edge(&graph, &dense_mixed_partitions, &sparse_mixed_partitions, indel_rate, bl)?;

    if metrics.second_derivative < 0.0 {
      let implied_step = (metrics.derivative / metrics.second_derivative).abs();
      let tol = newton_tolerance_t(bl);
      assert!(
        implied_step < tol * 10.0,
        "{method:?}: implied Newton step ({implied_step}) exceeds 10x tolerance ({tol}), \
         dl={}, d2l={}, t*={bl}",
        metrics.derivative, metrics.second_derivative
      );
    }

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::k1(1)]
  #[case::k2(2)]
  #[case::k4(4)]
  #[trace]
  fn test_optimize_method_cross_method_lh_agreement(#[case] n_indels: usize) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_brent_names = nwk_parsed.names();
    let graph_brent = nwk_parsed.graph;
    let mut bl_brent = nwk_parsed.branch_lengths;
    let (dense_partitions_brent, sparse_partitions_brent, rate_brent) = setup_with_indels(&graph_brent, &graph_brent_names, &mut bl_brent, n_indels)?;
    let total_length_brent = total_sequence_length(&dense_partitions_brent, &sparse_partitions_brent);
    let contributions_brent = gather_edge_contributions(&graph_brent, &dense_partitions_brent, &sparse_partitions_brent)?;
    let indel_counts_brent = gather_edge_indel_counts(&graph_brent, &dense_partitions_brent, &sparse_partitions_brent);
    run_optimize_mixed(&graph_brent, total_length_brent, &contributions_brent, &indel_counts_brent, BranchOptMethod::Brent, &mut bl_brent)?;
    let bl_brent = first_edge_bl(&graph_brent, &bl_brent);
    let lh_brent = eval_combined_first_edge(&graph_brent, &dense_partitions_brent, &sparse_partitions_brent, rate_brent, bl_brent)?;

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_sqrt_names = nwk_parsed.names();
    let graph_sqrt = nwk_parsed.graph;
    let mut bl_sqrt = nwk_parsed.branch_lengths;
    let (dense_partitions_sqrt, sparse_partitions_sqrt, rate_sqrt) = setup_with_indels(&graph_sqrt, &graph_sqrt_names, &mut bl_sqrt, n_indels)?;
    let total_length_sqrt = total_sequence_length(&dense_partitions_sqrt, &sparse_partitions_sqrt);
    let contributions_sqrt = gather_edge_contributions(&graph_sqrt, &dense_partitions_sqrt, &sparse_partitions_sqrt)?;
    let indel_counts_sqrt = gather_edge_indel_counts(&graph_sqrt, &dense_partitions_sqrt, &sparse_partitions_sqrt);
    run_optimize_mixed(&graph_sqrt, total_length_sqrt, &contributions_sqrt, &indel_counts_sqrt, BranchOptMethod::NewtonSqrt, &mut bl_sqrt)?;
    let bl_sqrt = first_edge_bl(&graph_sqrt, &bl_sqrt);
    let lh_sqrt = eval_combined_first_edge(&graph_sqrt, &dense_partitions_sqrt, &sparse_partitions_sqrt, rate_sqrt, bl_sqrt)?;

    let lh_diff = (lh_brent - lh_sqrt).abs();
    assert!(
      lh_diff < 1e-3,
      "Brent lh ({lh_brent}) and NewtonSqrt lh ({lh_sqrt}) differ by {lh_diff}, \
       bl_brent={bl_brent}, bl_sqrt={bl_sqrt}"
    );

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::k1(1)]
  #[case::k2(2)]
  #[case::k4(4)]
  #[trace]
  fn test_optimize_method_cross_method_lh_agreement_newton_log(#[case] n_indels: usize) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_brent_names = nwk_parsed.names();
    let graph_brent = nwk_parsed.graph;
    let mut bl_brent = nwk_parsed.branch_lengths;
    let (dense_partitions_brent, sparse_partitions_brent, rate_brent) = setup_with_indels(&graph_brent, &graph_brent_names, &mut bl_brent, n_indels)?;
    let total_length_brent = total_sequence_length(&dense_partitions_brent, &sparse_partitions_brent);
    let contributions_brent = gather_edge_contributions(&graph_brent, &dense_partitions_brent, &sparse_partitions_brent)?;
    let indel_counts_brent = gather_edge_indel_counts(&graph_brent, &dense_partitions_brent, &sparse_partitions_brent);
    run_optimize_mixed(&graph_brent, total_length_brent, &contributions_brent, &indel_counts_brent, BranchOptMethod::Brent, &mut bl_brent)?;
    let bl_brent = first_edge_bl(&graph_brent, &bl_brent);
    let lh_brent = eval_combined_first_edge(&graph_brent, &dense_partitions_brent, &sparse_partitions_brent, rate_brent, bl_brent)?;

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_log_names = nwk_parsed.names();
    let graph_log = nwk_parsed.graph;
    let mut bl_log = nwk_parsed.branch_lengths;
    let (dense_partitions_log, sparse_partitions_log, rate_log) = setup_with_indels(&graph_log, &graph_log_names, &mut bl_log, n_indels)?;
    let total_length_log = total_sequence_length(&dense_partitions_log, &sparse_partitions_log);
    let contributions_log = gather_edge_contributions(&graph_log, &dense_partitions_log, &sparse_partitions_log)?;
    let indel_counts_log = gather_edge_indel_counts(&graph_log, &dense_partitions_log, &sparse_partitions_log);
    run_optimize_mixed(&graph_log, total_length_log, &contributions_log, &indel_counts_log, BranchOptMethod::NewtonLog, &mut bl_log)?;
    let bl_log = first_edge_bl(&graph_log, &bl_log);
    let lh_log = eval_combined_first_edge(&graph_log, &dense_partitions_log, &sparse_partitions_log, rate_log, bl_log)?;

    let lh_diff = (lh_brent - lh_log).abs();
    assert!(
      lh_diff < 1e-3,
      "Brent lh ({lh_brent}) and NewtonLog lh ({lh_log}) differ by {lh_diff}, \
       bl_brent={bl_brent}, bl_log={bl_log}"
    );

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::newton_k1(     BranchOptMethod::Newton,      1)]
  #[case::newton_k2(     BranchOptMethod::Newton,      2)]
  #[case::newton_k4(     BranchOptMethod::Newton,      4)]
  #[case::newton_sqrt_k1(BranchOptMethod::NewtonSqrt,  1)]
  #[case::newton_sqrt_k2(BranchOptMethod::NewtonSqrt,  2)]
  #[case::newton_sqrt_k4(BranchOptMethod::NewtonSqrt,  4)]
  #[case::newton_log_k1( BranchOptMethod::NewtonLog,   1)]
  #[case::newton_log_k2( BranchOptMethod::NewtonLog,   2)]
  #[case::newton_log_k4( BranchOptMethod::NewtonLog,   4)]
  #[case::brent_k1(      BranchOptMethod::Brent,       1)]
  #[case::brent_k2(      BranchOptMethod::Brent,       2)]
  #[case::brent_k4(      BranchOptMethod::Brent,       4)]
  #[case::brent_log_k1(  BranchOptMethod::BrentLog,    1)]
  #[case::brent_log_k2(  BranchOptMethod::BrentLog,    2)]
  #[case::brent_log_k4(  BranchOptMethod::BrentLog,    4)]
  #[trace]
  fn test_optimize_method_cross_method_all_six_lh_agreement(
    #[case] method: BranchOptMethod,
    #[case] n_indels: usize,
  ) -> Result<(), Report> {
    let lh_ref = {
      let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
      let graph_ref_names = nwk_parsed.names();
      let graph_ref = nwk_parsed.graph;
      let mut bl_ref = nwk_parsed.branch_lengths;
      let (dense_partitions_ref, sparse_partitions_ref, rate_ref) = setup_with_indels(&graph_ref, &graph_ref_names, &mut bl_ref, n_indels)?;
      let total_length_ref = total_sequence_length(&dense_partitions_ref, &sparse_partitions_ref);
      let contributions_ref = gather_edge_contributions(&graph_ref, &dense_partitions_ref, &sparse_partitions_ref)?;
      let indel_counts_ref = gather_edge_indel_counts(&graph_ref, &dense_partitions_ref, &sparse_partitions_ref);
      run_optimize_mixed(&graph_ref, total_length_ref, &contributions_ref, &indel_counts_ref, BranchOptMethod::BrentSqrt, &mut bl_ref)?;
      let bl_ref = first_edge_bl(&graph_ref, &bl_ref);
      eval_combined_first_edge(&graph_ref, &dense_partitions_ref, &sparse_partitions_ref, rate_ref, bl_ref)?
    };

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let (dense_partitions, sparse_partitions, rate) = setup_with_indels(&graph, &names, &mut branch_lengths, n_indels)?;
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_partitions, &sparse_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);
    run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
    let bl = first_edge_bl(&graph, &branch_lengths);
    let lh = eval_combined_first_edge(&graph, &dense_partitions, &sparse_partitions, rate, bl)?;

    let diff = (lh - lh_ref).abs();
    assert!(
      diff < 1e-3,
      "{method:?} (n_indels={n_indels}) lh ({lh}) differs from BrentSqrt lh ({lh_ref}) by {diff} > 1e-3, bl={bl}"
    );

    Ok(())
  }

  #[test]
  fn test_optimize_method_newton_log_improves_over_newton() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_newton_names = nwk_parsed.names();
    let graph_newton = nwk_parsed.graph;
    let mut bl_newton = nwk_parsed.branch_lengths;
    let (dense_partitions_newton, sparse_partitions_newton, rate_newton) =
      setup_with_indels(&graph_newton, &graph_newton_names, &mut bl_newton, 4)?;
    let total_length_newton = total_sequence_length(&dense_partitions_newton, &sparse_partitions_newton);
    let contributions_newton =
      gather_edge_contributions(&graph_newton, &dense_partitions_newton, &sparse_partitions_newton)?;
    let indel_counts_newton =
      gather_edge_indel_counts(&graph_newton, &dense_partitions_newton, &sparse_partitions_newton);
    run_optimize_mixed(
      &graph_newton,
      total_length_newton,
      &contributions_newton,
      &indel_counts_newton,
      BranchOptMethod::Newton,
      &mut bl_newton,
    )?;
    let bl_newton = first_edge_bl(&graph_newton, &bl_newton);
    let lh_newton = eval_combined_first_edge(
      &graph_newton,
      &dense_partitions_newton,
      &sparse_partitions_newton,
      rate_newton,
      bl_newton,
    )?;

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_log_names = nwk_parsed.names();
    let graph_log = nwk_parsed.graph;
    let mut bl_log = nwk_parsed.branch_lengths;
    let (dense_partitions_log, sparse_partitions_log, rate_log) =
      setup_with_indels(&graph_log, &graph_log_names, &mut bl_log, 4)?;
    let total_length_log = total_sequence_length(&dense_partitions_log, &sparse_partitions_log);
    let contributions_log = gather_edge_contributions(&graph_log, &dense_partitions_log, &sparse_partitions_log)?;
    let indel_counts_log = gather_edge_indel_counts(&graph_log, &dense_partitions_log, &sparse_partitions_log);
    run_optimize_mixed(
      &graph_log,
      total_length_log,
      &contributions_log,
      &indel_counts_log,
      BranchOptMethod::NewtonLog,
      &mut bl_log,
    )?;
    let bl_log = first_edge_bl(&graph_log, &bl_log);
    let lh_log = eval_combined_first_edge(
      &graph_log,
      &dense_partitions_log,
      &sparse_partitions_log,
      rate_log,
      bl_log,
    )?;

    assert!(bl_newton > 0.0 && bl_newton.is_finite());
    assert!(bl_log > 0.0 && bl_log.is_finite());

    assert!(
      lh_log >= lh_newton - 1e-10,
      "NewtonLog lh ({lh_log}) should be >= Newton lh ({lh_newton}), \
       bl_log={bl_log}, bl_newton={bl_newton}"
    );

    Ok(())
  }

  #[test]
  fn test_optimize_method_newton_sqrt_improves_over_newton() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_newton_names = nwk_parsed.names();
    let graph_newton = nwk_parsed.graph;
    let mut bl_newton = nwk_parsed.branch_lengths;
    let (dense_partitions_newton, sparse_partitions_newton, rate_newton) =
      setup_with_indels(&graph_newton, &graph_newton_names, &mut bl_newton, 4)?;
    let total_length_newton = total_sequence_length(&dense_partitions_newton, &sparse_partitions_newton);
    let contributions_newton =
      gather_edge_contributions(&graph_newton, &dense_partitions_newton, &sparse_partitions_newton)?;
    let indel_counts_newton =
      gather_edge_indel_counts(&graph_newton, &dense_partitions_newton, &sparse_partitions_newton);
    run_optimize_mixed(
      &graph_newton,
      total_length_newton,
      &contributions_newton,
      &indel_counts_newton,
      BranchOptMethod::Newton,
      &mut bl_newton,
    )?;
    let bl_newton = first_edge_bl(&graph_newton, &bl_newton);
    let lh_newton = eval_combined_first_edge(
      &graph_newton,
      &dense_partitions_newton,
      &sparse_partitions_newton,
      rate_newton,
      bl_newton,
    )?;

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_sqrt_names = nwk_parsed.names();
    let graph_sqrt = nwk_parsed.graph;
    let mut bl_sqrt = nwk_parsed.branch_lengths;
    let (dense_partitions_sqrt, sparse_partitions_sqrt, rate_sqrt) =
      setup_with_indels(&graph_sqrt, &graph_sqrt_names, &mut bl_sqrt, 4)?;
    let total_length_sqrt = total_sequence_length(&dense_partitions_sqrt, &sparse_partitions_sqrt);
    let contributions_sqrt = gather_edge_contributions(&graph_sqrt, &dense_partitions_sqrt, &sparse_partitions_sqrt)?;
    let indel_counts_sqrt = gather_edge_indel_counts(&graph_sqrt, &dense_partitions_sqrt, &sparse_partitions_sqrt);
    run_optimize_mixed(
      &graph_sqrt,
      total_length_sqrt,
      &contributions_sqrt,
      &indel_counts_sqrt,
      BranchOptMethod::NewtonSqrt,
      &mut bl_sqrt,
    )?;
    let bl_sqrt = first_edge_bl(&graph_sqrt, &bl_sqrt);
    let lh_sqrt = eval_combined_first_edge(
      &graph_sqrt,
      &dense_partitions_sqrt,
      &sparse_partitions_sqrt,
      rate_sqrt,
      bl_sqrt,
    )?;

    assert!(bl_newton > 0.0 && bl_newton.is_finite());
    assert!(bl_sqrt > 0.0 && bl_sqrt.is_finite());

    assert!(
      lh_sqrt >= lh_newton - 1e-10,
      "NewtonSqrt lh ({lh_sqrt}) should be >= Newton lh ({lh_newton}), \
       bl_sqrt={bl_sqrt}, bl_newton={bl_newton}"
    );

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::k1(1)]
  #[case::k2(2)]
  #[case::k4(4)]
  #[trace]
  fn test_optimize_method_newton_cross_conditioning_ordering(#[case] n_indels: usize) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_newton_names = nwk_parsed.names();
    let graph_newton = nwk_parsed.graph;
    let mut bl_newton = nwk_parsed.branch_lengths;
    let (dense_partitions_newton, sparse_partitions_newton, rate_newton) = setup_with_indels(&graph_newton, &graph_newton_names, &mut bl_newton, n_indels)?;
    let total_length_newton = total_sequence_length(&dense_partitions_newton, &sparse_partitions_newton);
    let contributions_newton = gather_edge_contributions(&graph_newton, &dense_partitions_newton, &sparse_partitions_newton)?;
    let indel_counts_newton = gather_edge_indel_counts(&graph_newton, &dense_partitions_newton, &sparse_partitions_newton);
    run_optimize_mixed(&graph_newton, total_length_newton, &contributions_newton, &indel_counts_newton, BranchOptMethod::Newton, &mut bl_newton)?;
    let bl_newton = first_edge_bl(&graph_newton, &bl_newton);
    let lh_newton = eval_combined_first_edge(&graph_newton, &dense_partitions_newton, &sparse_partitions_newton, rate_newton, bl_newton)?;

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_sqrt_names = nwk_parsed.names();
    let graph_sqrt = nwk_parsed.graph;
    let mut bl_sqrt = nwk_parsed.branch_lengths;
    let (dense_partitions_sqrt, sparse_partitions_sqrt, rate_sqrt) = setup_with_indels(&graph_sqrt, &graph_sqrt_names, &mut bl_sqrt, n_indels)?;
    let total_length_sqrt = total_sequence_length(&dense_partitions_sqrt, &sparse_partitions_sqrt);
    let contributions_sqrt = gather_edge_contributions(&graph_sqrt, &dense_partitions_sqrt, &sparse_partitions_sqrt)?;
    let indel_counts_sqrt = gather_edge_indel_counts(&graph_sqrt, &dense_partitions_sqrt, &sparse_partitions_sqrt);
    run_optimize_mixed(&graph_sqrt, total_length_sqrt, &contributions_sqrt, &indel_counts_sqrt, BranchOptMethod::NewtonSqrt, &mut bl_sqrt)?;
    let bl_sqrt = first_edge_bl(&graph_sqrt, &bl_sqrt);
    let lh_sqrt = eval_combined_first_edge(&graph_sqrt, &dense_partitions_sqrt, &sparse_partitions_sqrt, rate_sqrt, bl_sqrt)?;

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_log_names = nwk_parsed.names();
    let graph_log = nwk_parsed.graph;
    let mut bl_log = nwk_parsed.branch_lengths;
    let (dense_partitions_log, sparse_partitions_log, rate_log) = setup_with_indels(&graph_log, &graph_log_names, &mut bl_log, n_indels)?;
    let total_length_log = total_sequence_length(&dense_partitions_log, &sparse_partitions_log);
    let contributions_log = gather_edge_contributions(&graph_log, &dense_partitions_log, &sparse_partitions_log)?;
    let indel_counts_log = gather_edge_indel_counts(&graph_log, &dense_partitions_log, &sparse_partitions_log);
    run_optimize_mixed(&graph_log, total_length_log, &contributions_log, &indel_counts_log, BranchOptMethod::NewtonLog, &mut bl_log)?;
    let bl_log = first_edge_bl(&graph_log, &bl_log);
    let lh_log = eval_combined_first_edge(&graph_log, &dense_partitions_log, &sparse_partitions_log, rate_log, bl_log)?;

    let tol = 1e-10;
    assert!(
      lh_sqrt >= lh_newton - tol,
      "NewtonSqrt lh ({lh_sqrt}) should be >= Newton lh ({lh_newton}), \
       bl_sqrt={bl_sqrt}, bl_newton={bl_newton}, k={n_indels}"
    );
    assert!(
      lh_log >= lh_sqrt - tol,
      "NewtonLog lh ({lh_log}) should be >= NewtonSqrt lh ({lh_sqrt}), \
       bl_log={bl_log}, bl_sqrt={bl_sqrt}, k={n_indels}"
    );

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::brent_k1(     BranchOptMethod::Brent,     1)]
  #[case::brent_k2(     BranchOptMethod::Brent,     2)]
  #[case::brent_k4(     BranchOptMethod::Brent,     4)]
  #[case::brent_sqrt_k1(BranchOptMethod::BrentSqrt, 1)]
  #[case::brent_sqrt_k2(BranchOptMethod::BrentSqrt, 2)]
  #[case::brent_sqrt_k4(BranchOptMethod::BrentSqrt, 4)]
  #[case::brent_log_k1( BranchOptMethod::BrentLog,  1)]
  #[case::brent_log_k2( BranchOptMethod::BrentLog,  2)]
  #[case::brent_log_k4( BranchOptMethod::BrentLog,  4)]
  #[trace]
  fn test_optimize_method_brent_positive_with_indels(
    #[case] method: BranchOptMethod,
    #[case] n_indels: usize,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (dense_mixed_partitions, sparse_mixed_partitions, _) = setup_with_indels(&graph, &names, &mut branch_lengths, n_indels)?;
    let total_length = total_sequence_length(&dense_mixed_partitions, &sparse_mixed_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_mixed_partitions, &sparse_mixed_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_mixed_partitions, &sparse_mixed_partitions);

    run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;

    let bl = first_edge_bl(&graph, &branch_lengths);
    assert!(bl > 0.0, "{method:?} BL with {n_indels} indels must be positive, got {bl}");
    assert!(bl.is_finite(), "{method:?} BL must be finite, got {bl}");
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::k1(1)]
  #[case::k2(2)]
  #[case::k4(4)]
  #[trace]
  fn test_optimize_method_newton_log_positive_with_indels(#[case] n_indels: usize) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (dense_mixed_partitions, sparse_mixed_partitions, _) = setup_with_indels(&graph, &names, &mut branch_lengths, n_indels)?;
    let total_length = total_sequence_length(&dense_mixed_partitions, &sparse_mixed_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_mixed_partitions, &sparse_mixed_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_mixed_partitions, &sparse_mixed_partitions);

    run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, BranchOptMethod::NewtonLog, &mut branch_lengths)?;

    let bl = first_edge_bl(&graph, &branch_lengths);
    assert!(bl > 0.0, "NewtonLog BL with {n_indels} indels must be positive, got {bl}");
    assert!(bl.is_finite(), "NewtonLog BL must be finite, got {bl}");
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::k1(1)]
  #[case::k2(2)]
  #[case::k4(4)]
  #[trace]
  fn test_optimize_method_newton_sqrt_positive_with_indels(#[case] n_indels: usize) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (dense_mixed_partitions, sparse_mixed_partitions, _) = setup_with_indels(&graph, &names, &mut branch_lengths, n_indels)?;
    let total_length = total_sequence_length(&dense_mixed_partitions, &sparse_mixed_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_mixed_partitions, &sparse_mixed_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_mixed_partitions, &sparse_mixed_partitions);

    run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, BranchOptMethod::NewtonSqrt, &mut branch_lengths)?;

    let bl = first_edge_bl(&graph, &branch_lengths);
    assert!(bl > 0.0, "NewtonSqrt BL with {n_indels} indels must be positive, got {bl}");
    assert!(bl.is_finite(), "NewtonSqrt BL must be finite, got {bl}");
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::k1(1)]
  #[case::k2(2)]
  #[case::k4(4)]
  #[trace]
  fn test_optimize_method_newton_positive_with_indels(#[case] n_indels: usize) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (dense_mixed_partitions, sparse_mixed_partitions, _) = setup_with_indels(&graph, &names, &mut branch_lengths, n_indels)?;
    let total_length = total_sequence_length(&dense_mixed_partitions, &sparse_mixed_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_mixed_partitions, &sparse_mixed_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_mixed_partitions, &sparse_mixed_partitions);

    run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, BranchOptMethod::Newton, &mut branch_lengths)?;

    let bl = first_edge_bl(&graph, &branch_lengths);
    assert!(bl > 0.0, "Newton BL with {n_indels} indels must be positive, got {bl}");
    assert!(bl.is_finite(), "Newton BL must be finite, got {bl}");
    Ok(())
  }

  #[test]
  fn test_optimize_method_brent_cross_parameterization_lh_agreement() -> Result<(), Report> {
    let n_indels = 3;

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_t_names = nwk_parsed.names();
    let graph_t = nwk_parsed.graph;
    let mut bl_t = nwk_parsed.branch_lengths;
    let (dense_parts_t, sparse_parts_t, rate_t) = setup_with_indels(&graph_t, &graph_t_names, &mut bl_t, n_indels)?;
    let total_length_t = total_sequence_length(&dense_parts_t, &sparse_parts_t);
    let contributions_t = gather_edge_contributions(&graph_t, &dense_parts_t, &sparse_parts_t)?;
    let indel_counts_t = gather_edge_indel_counts(&graph_t, &dense_parts_t, &sparse_parts_t);
    run_optimize_mixed(
      &graph_t,
      total_length_t,
      &contributions_t,
      &indel_counts_t,
      BranchOptMethod::Brent,
      &mut bl_t,
    )?;
    let bl_t = first_edge_bl(&graph_t, &bl_t);
    let lh_t = eval_combined_first_edge(&graph_t, &dense_parts_t, &sparse_parts_t, rate_t, bl_t)?;

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_sqrt_names = nwk_parsed.names();
    let graph_sqrt = nwk_parsed.graph;
    let mut bl_sqrt = nwk_parsed.branch_lengths;
    let (dense_parts_sqrt, sparse_parts_sqrt, rate_sqrt) =
      setup_with_indels(&graph_sqrt, &graph_sqrt_names, &mut bl_sqrt, n_indels)?;
    let total_length_sqrt = total_sequence_length(&dense_parts_sqrt, &sparse_parts_sqrt);
    let contributions_sqrt = gather_edge_contributions(&graph_sqrt, &dense_parts_sqrt, &sparse_parts_sqrt)?;
    let indel_counts_sqrt = gather_edge_indel_counts(&graph_sqrt, &dense_parts_sqrt, &sparse_parts_sqrt);
    run_optimize_mixed(
      &graph_sqrt,
      total_length_sqrt,
      &contributions_sqrt,
      &indel_counts_sqrt,
      BranchOptMethod::BrentSqrt,
      &mut bl_sqrt,
    )?;
    let bl_sqrt = first_edge_bl(&graph_sqrt, &bl_sqrt);
    let lh_sqrt = eval_combined_first_edge(&graph_sqrt, &dense_parts_sqrt, &sparse_parts_sqrt, rate_sqrt, bl_sqrt)?;

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_log_names = nwk_parsed.names();
    let graph_log = nwk_parsed.graph;
    let mut bl_log = nwk_parsed.branch_lengths;
    let (dense_parts_log, sparse_parts_log, rate_log) =
      setup_with_indels(&graph_log, &graph_log_names, &mut bl_log, n_indels)?;
    let total_length_log = total_sequence_length(&dense_parts_log, &sparse_parts_log);
    let contributions_log = gather_edge_contributions(&graph_log, &dense_parts_log, &sparse_parts_log)?;
    let indel_counts_log = gather_edge_indel_counts(&graph_log, &dense_parts_log, &sparse_parts_log);
    run_optimize_mixed(
      &graph_log,
      total_length_log,
      &contributions_log,
      &indel_counts_log,
      BranchOptMethod::BrentLog,
      &mut bl_log,
    )?;
    let bl_log = first_edge_bl(&graph_log, &bl_log);
    let lh_log = eval_combined_first_edge(&graph_log, &dense_parts_log, &sparse_parts_log, rate_log, bl_log)?;

    let diff_sqrt = (lh_t - lh_sqrt).abs();
    let diff_log = (lh_t - lh_log).abs();
    assert!(
      diff_sqrt < 1e-3,
      "Brent-t lh ({lh_t}) vs Brent-sqrt lh ({lh_sqrt}) differ by {diff_sqrt}"
    );
    assert!(
      diff_log < 1e-3,
      "Brent-t lh ({lh_t}) vs Brent-log lh ({lh_log}) differ by {diff_log}"
    );

    Ok(())
  }

  mod generators {
    use proptest::prelude::*;
    pub(crate) fn gen_s() -> impl Strategy<Value = f64> {
      1e-6_f64..1e3_f64
    }
    pub(crate) fn gen_t() -> impl Strategy<Value = f64> {
      1e-10_f64..1e3_f64
    }
    pub(crate) fn gen_dl_dt() -> impl Strategy<Value = f64> {
      -1e6_f64..1e6_f64
    }
    pub(crate) fn gen_d2l_dt2() -> impl Strategy<Value = f64> {
      -1e8_f64..1e8_f64
    }
    pub(crate) fn gen_scalar() -> impl Strategy<Value = f64> {
      prop_oneof![-1e3_f64..-1e-3_f64, 1e-3_f64..1e3_f64]
    }
  }

  pub(crate) mod helpers {
    use super::*;

    pub(crate) fn setup_with_indels(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
      n_indels: usize,
    ) -> Result<(Vec<DenseReconstruction>, Vec<SparseReconstruction>, f64), Report> {
      let aln = simple_alignment()?;
      let (mut dense_partitions, mut sparse_partitions) = setup_partitions(graph, names, &aln, branch_lengths)?;

      let first_edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
      let indels: Vec<InDel> = (0..n_indels)
        .map(|i| InDel::del((i * 3, i * 3 + 3), Seq::try_from_str("ACG").unwrap()).unwrap())
        .collect();

      for p in &mut dense_partitions {
        p.edges.estimates.get_mut(&first_edge_key).unwrap().indels = indels.clone();
      }
      for p in &mut sparse_partitions {
        p.partition.obs_edges.get_mut(&first_edge_key).unwrap().indels = indels.clone();
      }

      let indel_rate = {
        let indel_counts = gather_edge_indel_counts(graph, &dense_partitions, &sparse_partitions);
        estimate_indel_rate(graph, &indel_counts, branch_lengths)
      };

      Ok((dense_partitions, sparse_partitions, indel_rate))
    }

    pub(crate) fn eval_combined_first_edge(
      graph: &Graph,
      dense: &[DenseReconstruction],
      sparse: &[SparseReconstruction],
      indel_rate: f64,
      t: f64,
    ) -> Result<f64, Report> {
      let edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
      let contributions = gather_edge_contributions(graph, dense, sparse)?;
      let indel_counts = gather_edge_indel_counts(graph, dense, sparse);
      let edge_contributions = &contributions[&edge_key];
      let indel_count: usize = indel_counts[&edge_key];

      let sub_lh = evaluate_mixed_log_lh_only(edge_contributions, t)
        .expect("valid branch length")
        .value();
      let indel_lh = poisson_indel_log_lh(indel_count, indel_rate, t)
        .expect("valid Poisson parameters")
        .log_lh
        .value();
      Ok(sub_lh + indel_lh)
    }

    pub(crate) fn eval_metrics_first_edge(
      graph: &Graph,
      dense: &[DenseReconstruction],
      sparse: &[SparseReconstruction],
      indel_rate: f64,
      t: f64,
    ) -> Result<OptimizationMetrics, Report> {
      let edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
      let contributions = gather_edge_contributions(graph, dense, sparse)?;
      let indel_counts = gather_edge_indel_counts(graph, dense, sparse);
      let edge_contributions = &contributions[&edge_key];
      let indel_count: usize = indel_counts[&edge_key];

      let mut metrics = evaluate_mixed(edge_contributions, t).expect("valid branch length");
      metrics.add(&poisson_indel_log_lh(indel_count, indel_rate, t).expect("valid Poisson parameters"));
      Ok(metrics)
    }

    pub(crate) fn first_edge_bl(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> f64 {
      branch_lengths[&graph.get_edges().collect::<Vec<_>>()[0].key()].unwrap()
    }

    proptest! {
      #[test]
      fn test_prop_optimize_method_chain_rule_sqrt_formula(
        s in generators::gen_s(),
        dl_dt in generators::gen_dl_dt(),
        d2l_dt2 in generators::gen_d2l_dt2(),
      ) {
        let (dl_ds, d2l_ds2) = chain_rule_sqrt(s, dl_dt, d2l_dt2);
        let expected_dl_ds = 2.0 * s * dl_dt;
        let expected_d2l_ds2 = 4.0 * s * s * d2l_dt2 + 2.0 * dl_dt;
        let dl_tol = 1e-9 * expected_dl_ds.abs().max(1.0);
        let d2l_tol = 1e-9 * expected_d2l_ds2.abs().max(1.0);
        prop_assert!(
          (dl_ds - expected_dl_ds).abs() <= dl_tol,
          "dl_ds: got {dl_ds}, expected {expected_dl_ds}"
        );
        prop_assert!(
          (d2l_ds2 - expected_d2l_ds2).abs() <= d2l_tol,
          "d2l_ds2: got {d2l_ds2}, expected {expected_d2l_ds2}"
        );
      }

      #[test]
      fn test_prop_optimize_method_chain_rule_sqrt_linear(
        s in generators::gen_s(),
        dl_dt in generators::gen_dl_dt(),
        d2l_dt2 in generators::gen_d2l_dt2(),
        k in generators::gen_scalar(),
      ) {
        let (dl_ds, d2l_ds2) = chain_rule_sqrt(s, dl_dt, d2l_dt2);
        let (dl_ds_k, d2l_ds2_k) = chain_rule_sqrt(s, k * dl_dt, k * d2l_dt2);
        let dl_tol = 1e-9 * (k * dl_ds).abs().max(1.0);
        let d2l_tol = 1e-9 * (k * d2l_ds2).abs().max(1.0);
        prop_assert!(
          (dl_ds_k - k * dl_ds).abs() <= dl_tol,
          "linearity violated for dl_ds: {dl_ds_k} vs k*{dl_ds}"
        );
        prop_assert!(
          (d2l_ds2_k - k * d2l_ds2).abs() <= d2l_tol,
          "linearity violated for d2l_ds2: {d2l_ds2_k} vs k*{d2l_ds2}"
        );
      }

      #[test]
      fn test_prop_optimize_method_chain_rule_log_formula(
        t in generators::gen_t(),
        dl_dt in generators::gen_dl_dt(),
        d2l_dt2 in generators::gen_d2l_dt2(),
      ) {
        let (dl_du, d2l_du2) = chain_rule_log(t, dl_dt, d2l_dt2);
        let expected_dl_du = t * dl_dt;
        let expected_d2l_du2 = t * t * d2l_dt2 + t * dl_dt;
        let dl_tol = 1e-9 * expected_dl_du.abs().max(1.0);
        let d2l_tol = 1e-9 * expected_d2l_du2.abs().max(1.0);
        prop_assert!(
          (dl_du - expected_dl_du).abs() <= dl_tol,
          "dl_du: got {dl_du}, expected {expected_dl_du}"
        );
        prop_assert!(
          (d2l_du2 - expected_d2l_du2).abs() <= d2l_tol,
          "d2l_du2: got {d2l_du2}, expected {expected_d2l_du2}"
        );
      }

      #[test]
      fn test_prop_optimize_method_chain_rule_log_linear(
        t in generators::gen_t(),
        dl_dt in generators::gen_dl_dt(),
        d2l_dt2 in generators::gen_d2l_dt2(),
        k in generators::gen_scalar(),
      ) {
        let (dl_du, d2l_du2) = chain_rule_log(t, dl_dt, d2l_dt2);
        let (dl_du_k, d2l_du2_k) = chain_rule_log(t, k * dl_dt, k * d2l_dt2);
        let dl_tol = 1e-9 * (k * dl_du).abs().max(1.0);
        let d2l_tol = 1e-9 * (k * d2l_du2).abs().max(1.0);
        prop_assert!(
          (dl_du_k - k * dl_du).abs() <= dl_tol,
          "linearity violated for dl_du: {dl_du_k} vs k*{dl_du}"
        );
        prop_assert!(
          (d2l_du2_k - k * d2l_du2).abs() <= d2l_tol,
          "linearity violated for d2l_du2: {d2l_du2_k} vs k*{d2l_du2}"
        );
      }
    }
  }
}
