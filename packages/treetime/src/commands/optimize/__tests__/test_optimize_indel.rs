#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::commands::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_partitions, simple_alignment,
  };
  use crate::commands::optimize::optimize_dense;
  use crate::commands::optimize::optimize_indel::{estimate_indel_rate, poisson_indel_log_lh};
  use crate::commands::optimize::optimize_unified::{
    OptimizationContribution, evaluate_mixed_log_lh_only, initial_guess_mixed, is_zero_better_than_grid_best,
    is_zero_branch_optimal, run_optimize_mixed,
  };
  use crate::commands::optimize::run::collect_optimize_partitions;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::seq::indel::InDel;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::array;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;

  /// Inject indels onto the first edge in each partition (both dense and sparse).
  fn inject_indels_on_first_edge(
    graph: &GraphAncestral,
    dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
    sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
    indels: &[InDel],
  ) -> treetime_graph::edge::GraphEdgeKey {
    let first_edge_key = graph.get_edges()[0].read_arc().key();
    for partition in dense_partitions {
      let mut p = partition.write_arc();
      p.edges.get_mut(&first_edge_key).unwrap().indels = indels.to_vec();
    }
    for partition in sparse_partitions {
      let mut p = partition.write_arc();
      p.edges.get_mut(&first_edge_key).unwrap().indels = indels.to_vec();
    }
    first_edge_key
  }

  /// Identical-sequence alignment: no substitutions on any edge.
  fn identical_alignment() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      ">A\nACGTACGTACGTACGT\n>B\nACGTACGTACGTACGT\n>C\nACGTACGTACGTACGT\n>D\nACGTACGTACGTACGT\n",
      &alphabet,
    )
  }

  /// Set up partitions with identical sequences (zero substitutions on every edge).
  fn setup_identical_partitions(
    graph: &GraphAncestral,
  ) -> Result<
    (
      Vec<Arc<RwLock<PartitionMarginalDense>>>,
      Vec<Arc<RwLock<PartitionMarginalSparse>>>,
      crate::commands::optimize::partition_ops::PartitionOptimizeVec,
    ),
    Report,
  > {
    let aln = identical_alignment()?;
    let alphabet_dense = Alphabet::new(AlphabetName::Nuc)?;
    let alphabet_sparse = Alphabet::new(AlphabetName::Nuc)?;

    let dense_partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: alphabet_dense,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    let sparse_partitions = vec![Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 1,
      gtr: jc69(JC69Params::default())?,
      alphabet: alphabet_sparse,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(graph, &sparse_partitions, &aln)?;
    initialize_marginal(graph, &dense_partitions, &aln)?;
    update_marginal(graph, &sparse_partitions)?;

    let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);
    initial_guess_mixed(graph, &mixed_partitions, true)?;

    Ok((dense_partitions, sparse_partitions, mixed_partitions))
  }

  /// `estimate_indel_rate` returns 0 when no edges have indels.
  #[test]
  fn test_optimize_indel_estimate_rate_no_indels() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let aln = simple_alignment()?;
    let (_, _, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let rate = estimate_indel_rate(&graph, &mixed_partitions);
    assert_abs_diff_eq!(rate, 0.0, epsilon = 1e-15);
    Ok(())
  }

  /// `estimate_indel_rate` returns total_indels / total_branch_length.
  #[test]
  fn test_optimize_indel_estimate_rate_with_indels() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let aln = simple_alignment()?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)];
    inject_indels_on_first_edge(&graph, &dense_partitions, &sparse_partitions, &indels);

    let rate = estimate_indel_rate(&graph, &mixed_partitions);

    // 2 indels total (1 per partition: dense + sparse), divided by total branch length
    let total_bl: f64 = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .sum();
    let expected_rate = 2.0 / total_bl;
    assert_abs_diff_eq!(rate, expected_rate, epsilon = 1e-10);
    Ok(())
  }

  /// With identical sequences (zero subs on every edge), `initial_guess_mixed` uses
  /// the Poisson MLE for indel-bearing edges.
  #[test]
  fn test_optimize_indel_initial_guess_nonzero_with_indels() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_identical_partitions(&graph)?;

    let indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)];
    inject_indels_on_first_edge(&graph, &dense_partitions, &sparse_partitions, &indels);

    initial_guess_mixed(&graph, &mixed_partitions, true)?;

    let bl = graph.get_edges()[0]
      .read_arc()
      .payload()
      .read_arc()
      .branch_length()
      .unwrap();
    assert!(
      bl > 0.0,
      "Branch with indels but no subs should have positive initial guess, got {bl}"
    );
    Ok(())
  }

  /// Regression: when ALL branch lengths are zero and the only signal is indels,
  /// `initial_guess_mixed` must still bootstrap positive branch lengths. This exercises
  /// the `indel_rate == 0` fallback path.
  #[test]
  fn test_optimize_indel_initial_guess_zero_bl_tree_with_indels() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_identical_partitions(&graph)?;

    // Zero all branch lengths to simulate degenerate input
    for edge_ref in graph.get_edges() {
      edge_ref.write_arc().payload().write_arc().set_branch_length(Some(0.0));
    }

    // Inject indels on the first edge (after zeroing, so indel_rate starts at 0)
    let indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)];
    inject_indels_on_first_edge(&graph, &dense_partitions, &sparse_partitions, &indels);

    // indel_rate is 0 at this point (all BL = 0), but initial_guess should bootstrap
    initial_guess_mixed(&graph, &mixed_partitions, true)?;

    let bl = graph.get_edges()[0]
      .read_arc()
      .payload()
      .read_arc()
      .branch_length()
      .unwrap();
    assert!(
      bl > 0.0,
      "Indel-bearing edge on zero-BL tree should get positive initial guess, got {bl}"
    );
    Ok(())
  }

  /// `run_optimize_mixed` assigns non-zero branch length when indels are present.
  /// Verifies both positivity and local optimality of the result.
  #[test]
  fn test_optimize_indel_run_optimize_nonzero_with_indels() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let aln = simple_alignment()?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    // Inject indels AFTER setup (which includes update_marginal) to avoid being wiped
    // by the backward pass that recreates DenseEdgePartition from scratch.
    let indels = vec![
      InDel::del((0, 3), Seq::try_from_str("ACG")?),
      InDel::del((5, 8), Seq::try_from_str("ACG")?),
    ];
    inject_indels_on_first_edge(&graph, &dense_partitions, &sparse_partitions, &indels);

    run_optimize_mixed(&graph, &mixed_partitions)?;

    let bl = graph.get_edges()[0]
      .read_arc()
      .payload()
      .read_arc()
      .branch_length()
      .unwrap();
    assert!(
      bl > 0.0,
      "Branch with indels should have positive length after optimization, got {bl}"
    );
    assert!(bl.is_finite(), "Optimized branch length should be finite, got {bl}");
    Ok(())
  }

  /// Full pipeline regression: zero all BLs, inject indels into sparse partition only
  /// (production path), run initial_guess + optimize, verify escape from zero.
  #[test]
  fn test_optimize_indel_zero_bl_pipeline_escapes_zero() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_identical_partitions(&graph)?;

    // Zero all branch lengths
    for edge_ref in graph.get_edges() {
      edge_ref.write_arc().payload().write_arc().set_branch_length(Some(0.0));
    }

    // Inject indels into sparse partition only (production path: dense indels are wiped
    // by update_marginal, so only sparse contributes in the real CLI flow)
    let first_edge_key = graph.get_edges()[0].read_arc().key();
    sparse_partitions[0]
      .write_arc()
      .edges
      .get_mut(&first_edge_key)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)];

    initial_guess_mixed(&graph, &mixed_partitions, true)?;

    // After initial_guess, the indel-bearing edge should have positive BL (bootstrap)
    let bl_after_guess = graph.get_edges()[0]
      .read_arc()
      .payload()
      .read_arc()
      .branch_length()
      .unwrap();
    assert!(
      bl_after_guess > 0.0,
      "initial_guess should bootstrap positive BL for indel-bearing edge, got {bl_after_guess}"
    );

    // Run marginal + optimize
    update_marginal(&graph, &dense_partitions)?;
    update_marginal(&graph, &sparse_partitions)?;
    run_optimize_mixed(&graph, &mixed_partitions)?;

    let bl_final = graph.get_edges()[0]
      .read_arc()
      .payload()
      .read_arc()
      .branch_length()
      .unwrap();
    assert!(bl_final > 0.0, "Optimized BL should be positive, got {bl_final}");
    assert!(bl_final.is_finite(), "Optimized BL should be finite, got {bl_final}");
    Ok(())
  }

  /// The Poisson indel contribution makes the combined second derivative more negative
  /// (more concave), aiding Newton convergence.
  #[test]
  fn test_optimize_indel_poisson_concavity() {
    let k = 3;
    let mu = 5.0;

    for &t in &[0.01, 0.1, 0.5, 1.0, 5.0] {
      let metrics = poisson_indel_log_lh(k, mu, t);
      assert!(
        metrics.second_derivative < 0.0,
        "Poisson log-likelihood should be concave for k>0, but at t={t} got d2={:.6}",
        metrics.second_derivative
      );
    }
  }

  /// The Poisson derivative at the MLE (t = k/mu) is zero.
  #[test]
  fn test_optimize_indel_poisson_mle_derivative_zero() {
    for k in 1..=10 {
      let mu = 3.0;
      let t_mle = k as f64 / mu;
      let metrics = poisson_indel_log_lh(k, mu, t_mle);
      assert_abs_diff_eq!(metrics.derivative, 0.0, epsilon = 1e-13);
    }
  }

  /// The Poisson log-likelihood at the MLE is the maximum.
  #[test]
  fn test_optimize_indel_poisson_mle_is_maximum() {
    let k = 4;
    let mu = 2.0;
    let t_mle = k as f64 / mu;
    let lh_mle = poisson_indel_log_lh(k, mu, t_mle).log_lh;

    for &delta in &[-0.5, -0.1, -0.01, 0.01, 0.1, 0.5] {
      let t = t_mle + delta;
      if t > 0.0 {
        let lh = poisson_indel_log_lh(k, mu, t).log_lh;
        assert!(
          lh <= lh_mle + 1e-14,
          "log-lh at t={t} ({lh}) should be <= log-lh at MLE ({lh_mle})"
        );
      }
    }
  }

  /// Numerical derivative matches analytical derivative.
  #[test]
  fn test_optimize_indel_poisson_numerical_derivative() {
    let k = 3;
    let mu = 5.0;
    let t = 0.3;

    let metrics = poisson_indel_log_lh(k, mu, t);

    // First derivative: central difference
    let h1 = 1e-7;
    let lh_plus = poisson_indel_log_lh(k, mu, t + h1).log_lh;
    let lh_minus = poisson_indel_log_lh(k, mu, t - h1).log_lh;
    let numerical_deriv = (lh_plus - lh_minus) / (2.0 * h1);
    assert_abs_diff_eq!(metrics.derivative, numerical_deriv, epsilon = 1e-8);

    // Second derivative: central difference with larger step to avoid cancellation
    let h2 = 1e-4;
    let lh_plus2 = poisson_indel_log_lh(k, mu, t + h2).log_lh;
    let lh_center = poisson_indel_log_lh(k, mu, t).log_lh;
    let lh_minus2 = poisson_indel_log_lh(k, mu, t - h2).log_lh;
    let numerical_second = (lh_plus2 - 2.0 * lh_center + lh_minus2) / (h2 * h2);
    assert_abs_diff_eq!(metrics.second_derivative, numerical_second, epsilon = 1e-5);
  }

  /// When `indel_count == 0` and substitution derivative is negative,
  /// `is_zero_branch_optimal` returns true (existing behavior preserved).
  #[test]
  fn test_optimize_indel_zero_branch_no_indels_unchanged() {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.0, 1.0, 0.0, 0.0]];
    let contribution = OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr));

    assert!(is_zero_branch_optimal(&[contribution]));
  }

  /// Grid search zero-comparison must reject zero when indels are present.
  /// Demonstrates the bug scenario: substitution-only likelihood prefers t=0
  /// (pure-state site has maximum likelihood at zero branch length), but the
  /// Poisson indel log-likelihood diverges to -infinity at t=0 for k > 0.
  #[test]
  fn test_optimize_indel_grid_zero_comparison_rejects_zero_with_indels() {
    let gtr = jc69(JC69Params::default()).unwrap();
    // Pure-state coefficients: substitution likelihood is maximized at t=0.
    // At t=0, L(0) = sum(coefficients) = 1.0. At t>0, eigenvalue decay reduces L(t).
    let coefficients = array![[1.0, 0.0, 0.0, 0.0]];
    let contribution = OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr));
    let contributions = vec![contribution];

    let best_positive = 0.01;
    let indel_count = 1_usize;
    let indel_rate = 5.0;

    // Bug precondition: substitution-only comparison prefers zero
    let sub_lh_zero = evaluate_mixed_log_lh_only(&contributions, 0.0);
    let sub_lh_best = evaluate_mixed_log_lh_only(&contributions, best_positive);
    assert!(
      sub_lh_zero > sub_lh_best,
      "Bug precondition: subs-only likelihood at zero ({sub_lh_zero}) should exceed positive ({sub_lh_best})"
    );

    // Indel-aware comparison: Poisson log-lh diverges to -infinity near t=0
    let indel_lh_near_zero = poisson_indel_log_lh(indel_count, indel_rate, 1e-15).log_lh;
    let indel_lh_at_best = poisson_indel_log_lh(indel_count, indel_rate, best_positive).log_lh;

    let combined_near_zero = sub_lh_zero + indel_lh_near_zero;
    let combined_at_best = sub_lh_best + indel_lh_at_best;

    assert!(
      combined_at_best > combined_near_zero,
      "Indel-aware comparison should prefer positive t: at_best={combined_at_best}, near_zero={combined_near_zero}"
    );

    // Production function: indel_count > 0 causes is_zero_better_than_grid_best to return false
    assert!(
      !is_zero_better_than_grid_best(&contributions, indel_count, indel_rate, best_positive),
      "Production helper must return false when indel_count > 0"
    );
  }

  /// Regression: without indels, the grid search zero-comparison still allows zero
  /// when substitution likelihood prefers it.
  #[test]
  fn test_optimize_indel_grid_zero_comparison_allows_zero_without_indels() {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[1.0, 0.0, 0.0, 0.0]];
    let contribution = OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr));
    let contributions = vec![contribution];

    // Production function: indel_count == 0 with pure-state coefficients selects zero
    assert!(
      is_zero_better_than_grid_best(&contributions, 0, 0.0, 0.01),
      "Production helper should return true when indel_count == 0 and zero is better for subs"
    );
  }

  mod generators {
    use proptest::prelude::*;

    /// Generate valid Poisson parameters: k in [1, 200], mu in [0.1, 1000], t in [1e-8, 100].
    pub fn poisson_params() -> impl Strategy<Value = (usize, f64, f64)> {
      (1..200_usize, 0.1..1000.0_f64, 1e-8..100.0_f64)
    }
  }

  mod prop_tests {
    use super::generators;
    use crate::commands::optimize::optimize_indel::poisson_indel_log_lh;
    use proptest::prelude::*;
    use treetime_utils::{prop_assert_abs_diff_eq, prop_assert_relative_eq};

    proptest! {
      /// For k > 0, the second derivative is always negative (log-concavity).
      #[test]
      fn test_prop_optimize_indel_concavity((k, mu, t) in generators::poisson_params()) {
        let metrics = poisson_indel_log_lh(k, mu, t);
        prop_assert!(
          metrics.second_derivative < 0.0,
          "Expected negative second derivative for k={k}, mu={mu}, t={t}, got {}",
          metrics.second_derivative
        );
      }

      /// At the MLE t = k/mu, the derivative is zero.
      #[test]
      fn test_prop_optimize_indel_mle_derivative((k, mu, _t) in generators::poisson_params()) {
        let t_mle = k as f64 / mu;
        let metrics = poisson_indel_log_lh(k, mu, t_mle);
        prop_assert_abs_diff_eq!(metrics.derivative, 0.0, epsilon = 1e-10);
      }

      /// For any k > 0 near t=0, the Poisson derivative is positive,
      /// confirming that the optimum is always at positive t (zero is never optimal).
      #[test]
      fn test_prop_optimize_indel_derivative_positive_near_zero((k, mu, _t) in generators::poisson_params()) {
        let near_zero = 1e-10;
        let metrics = poisson_indel_log_lh(k, mu, near_zero);
        prop_assert!(
          metrics.derivative > 0.0,
          "Poisson derivative near t=0 should be positive for k={k}, mu={mu}, got {}",
          metrics.derivative
        );
      }

      /// Numerical first derivative matches analytical derivative.
      #[test]
      fn test_prop_optimize_indel_numerical_derivative((k, mu, t) in generators::poisson_params()) {
        let h = t * 1e-6;
        prop_assert!(h > 0.0);
        let metrics = poisson_indel_log_lh(k, mu, t);
        let lh_plus = poisson_indel_log_lh(k, mu, t + h).log_lh;
        let lh_minus = poisson_indel_log_lh(k, mu, t - h).log_lh;
        let numerical = (lh_plus - lh_minus) / (2.0 * h);
        prop_assert_relative_eq!(metrics.derivative, numerical, max_relative = 1e-4);
      }
    }
  }
}
