#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{GtrModelName, JC69Params, get_gtr_by_name, jc69};
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::run_loop::find_zero_optimal_internal_edges;
  use crate::optimize::method_newton::{newton_inner, newton_sqrt_inner};
  use crate::optimize::dispatch::{initial_guess_mixed, run_optimize_mixed};
  use crate::optimize::likelihood::{evaluate_mixed, evaluate_mixed_log_lh_only};
  use crate::optimize::zero_boundary::{is_zero_branch_optimal, reconcile_zero_boundary};
  use crate::optimize::run_loop::collect_optimize_partitions;
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::partition::marginal_dense::PartitionMarginalDense;
  use crate::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::partition::optimization_contribution::OptimizationContribution;
  use crate::partition::optimize_dense;
  use crate::partition::traits::PartitionOptimizeVec;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use indoc::indoc;
  
  use ndarray::array;
  use parking_lot::RwLock;
  use rstest::rstest;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  /// 4-taxon tree with positive initial branch lengths on every edge.
  ///
  /// The optimizer therefore starts from $t > 0$ and each method must reach
  /// $t = 0$ to demonstrate the dispatch-level boundary check (not the
  /// initial-guess zero pass-through).
  const IDENTICAL_TREE_NEWICK: &str = "((A:0.1,B:0.1)AB:0.1,(C:0.1,D:0.1)CD:0.1)root:0.1;";

  /// Four identical sequences. Under any reversible model, the per-site
  /// likelihood of staying in the observed state is strictly decreasing in
  /// $t$, so the global ML optimum on every edge is $t = 0$.
  const IDENTICAL_ALIGNMENT: &str = indoc! {r#"
    >A
    ACGTACGTACGTACGT
    >B
    ACGTACGTACGTACGT
    >C
    ACGTACGTACGTACGT
    >D
    ACGTACGTACGTACGT
  "#};

  /// Set up dense + sparse partitions using the given GTR model with identical
  /// leaf sequences. `initial_guess_mixed` is called with `overwrite_valid =
  /// false` so the positive branch lengths from the newick string survive
  /// into `run_optimize_mixed`, forcing each method to perform real
  /// optimization work.
  fn setup_identical_partitions(
    graph: &GraphAncestral,
    model: GtrModelName,
  ) -> Result<
    (
      Vec<Arc<RwLock<PartitionMarginalDense>>>,
      Vec<Arc<RwLock<PartitionMarginalSparse>>>,
      PartitionOptimizeVec,
    ),
    Report,
  > {
    let aln = read_many_fasta_str(IDENTICAL_ALIGNMENT, &Alphabet::default())?;

    let dense_partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense::new(0, get_gtr_by_name(model)?, Alphabet::new(AlphabetName::Nuc)?, get_common_length(&aln)?)))];

    let fitch = create_fitch_partition(graph, 1, Alphabet::new(AlphabetName::Nuc)?, &aln)?;
    let sparse_partitions = vec![Arc::new(RwLock::new(
      fitch.into_marginal_sparse(get_gtr_by_name(model)?, graph)?,
    ))];
    initialize_marginal(graph, &dense_partitions, &aln)?;
    update_marginal(graph, &sparse_partitions)?;

    let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);
    initial_guess_mixed(graph, &mixed_partitions, false)?;

    Ok((dense_partitions, sparse_partitions, mixed_partitions))
  }

  /// All 6 methods must return exactly $t = 0$ on every edge for identical
  /// sequences under K80.
  ///
  /// K80 is not proven unimodal (Dinh and Matsen 2017, Corollary 3.1 does
  /// not apply), so the pre-dispatch `is_zero_branch_optimal` shortcut is
  /// bypassed. NewtonLog ($u = \ln t$) and the Brent variants cannot
  /// evaluate at $t = 0$ because their internal domain is strictly
  /// positive, so they would otherwise bias the optimum to a tiny positive
  /// value like $10^{-12}$. The post-dispatch boundary check in
  /// `run_optimize_mixed` must compare the optimizer's result against
  /// $t = 0$ and replace it when zero has a strictly higher log-likelihood.
  ///
  /// Downstream, `find_zero_optimal_internal_edges()` collapses internal
  /// edges with $t = 0$ during topology cleanup. Without this check the
  /// optimize loop's topology cleanup misses these edges for non-unimodal
  /// models.
  #[rustfmt::skip]
  #[rstest]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[case::brent(      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethod::BrentLog)]
  #[trace]
  fn test_dispatch_zero_boundary_k80_identical_sequences(
    #[case] method: BranchOptMethod,
  ) -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(IDENTICAL_TREE_NEWICK)?;
    let (_, _, mixed_partitions) = setup_identical_partitions(&graph, GtrModelName::K80)?;

    // Sanity check: every edge must start with a positive branch length so
    // that the optimizer has real work to do. If this fails, the test
    // collapses to a trivial pass-through of the initial guess.
    for (i, edge_ref) in graph.get_edges().iter().enumerate() {
      let bl = edge_ref.read_arc().payload().read_arc().branch_length().unwrap();
      assert!(
        bl > 0.0,
        "precondition: edge {i} must start with positive BL, got {bl}"
      );
    }

    run_optimize_mixed(&graph, &mixed_partitions, method)?;

    for (i, edge_ref) in graph.get_edges().iter().enumerate() {
      let bl = edge_ref.read_arc().payload().read_arc().branch_length().unwrap();
      assert!(
        bl == 0.0,
        "{method:?}: edge {i} branch length must be exactly 0 after optimization, got {bl}"
      );
    }

    Ok(())
  }

  /// Reachability for every non-unimodal nucleotide model the code
  /// classifies as multimodal: on identical sequences, `BrentSqrt` (the
  /// default method, which cannot evaluate exactly at $t = 0$) must still
  /// return zero on every edge after optimization. This exercises the
  /// post-dispatch reconciliation path across the full set of models for
  /// which the pre-dispatch `is_zero_branch_optimal` shortcut is bypassed,
  /// catching a bug that happens to be K80-only and does not manifest on
  /// the other non-unimodal models.
  ///
  /// Jtt92 is excluded because it operates on a 20-state amino acid
  /// alphabet; the test alignment is nucleotide.
  #[rustfmt::skip]
  #[rstest]
  #[case::k80(  GtrModelName::K80)]
  #[case::hky85(GtrModelName::HKY85)]
  #[case::t92(  GtrModelName::T92)]
  #[case::tn93( GtrModelName::TN93)]
  #[trace]
  fn test_dispatch_zero_boundary_non_unimodal_models_all_reach_zero(
    #[case] model: GtrModelName,
  ) -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(IDENTICAL_TREE_NEWICK)?;
    let (dense_partitions, _, mixed_partitions) = setup_identical_partitions(&graph, model)?;

    // Precondition: the model must be classified as non-unimodal so that
    // the pre-dispatch shortcut is bypassed and the post-dispatch
    // reconciliation is exercised.
    assert!(
      !dense_partitions[0].read_arc().data.gtr.unimodal_branch_likelihood,
      "precondition: {model:?} must be classified as non-unimodal"
    );

    run_optimize_mixed(&graph, &mixed_partitions, BranchOptMethod::BrentSqrt)?;

    for (i, edge_ref) in graph.get_edges().iter().enumerate() {
      let bl = edge_ref.read_arc().payload().read_arc().branch_length().unwrap();
      assert!(
        bl == 0.0,
        "{model:?}: edge {i} branch length must be exactly 0 after optimization, got {bl}"
      );
    }

    Ok(())
  }

  /// Pre-dispatch shortcut path: for unimodal models (JC69 in this test),
  /// `is_zero_branch_optimal` must return true on identical-sequence
  /// contributions, and the full `run_optimize_mixed` flow must also
  /// return zero through that shortcut. This proves that the shortcut
  /// itself reaches the correct answer, independently from the
  /// post-dispatch reconciliation path exercised by the non-unimodal
  /// test above.
  ///
  /// Verifying both paths on the same identical-sequences fixture also
  /// establishes that the fix is behaviorally transparent for unimodal
  /// models: a change in which path fires must not affect the result.
  #[test]
  fn test_dispatch_zero_boundary_jc69_pre_dispatch_shortcut_reaches_zero() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(IDENTICAL_TREE_NEWICK)?;
    let (dense_partitions, _, mixed_partitions) = setup_identical_partitions(&graph, GtrModelName::JC69)?;

    assert!(
      dense_partitions[0].read_arc().data.gtr.unimodal_branch_likelihood,
      "precondition: JC69 must be classified as unimodal"
    );

    // Build the contributions for the first edge and assert the
    // pre-dispatch derivative shortcut fires. Identical sequences produce
    // a strictly negative derivative at $t = 0$ for any unimodal model,
    // so the shortcut must return true.
    let first_edge_key = graph.get_edges()[0].read_arc().key();
    let contributions: Vec<OptimizationContribution> = mixed_partitions
      .iter()
      .map(|partition| partition.read_arc().create_edge_contribution(first_edge_key))
      .collect::<Result<_, _>>()?;
    assert!(
      is_zero_branch_optimal(&contributions),
      "precondition: JC69 identical-sequence contributions must trigger the pre-dispatch shortcut"
    );

    // Full dispatch also reaches zero through the shortcut path.
    run_optimize_mixed(&graph, &mixed_partitions, BranchOptMethod::BrentSqrt)?;
    for (i, edge_ref) in graph.get_edges().iter().enumerate() {
      let bl = edge_ref.read_arc().payload().read_arc().branch_length().unwrap();
      assert!(
        bl == 0.0,
        "JC69: edge {i} branch length must be exactly 0 after optimization, got {bl}"
      );
    }

    Ok(())
  }

  /// Build the Dinh and Matsen 2017 K80 $\kappa = 3$ counterexample as an
  /// `OptimizationContribution`. The eigenvalues and coefficients reproduce
  /// Section 5 (eq 5.1, eq 5.2) of the paper, also used in
  /// `test_is_zero_branch_optimal_k80_dinh_matsen_multimodal_counterexample`.
  /// The resulting surface has a positive local maximum near $t \approx 0.2$
  /// with $\ell(0.2) > \ell(0)$.
  fn make_dinh_matsen_k80_contribution() -> OptimizationContribution {
    let mut gtr = jc69(JC69Params::default()).unwrap();
    gtr.eigvals = array![-1.0, -0.5, -0.5, 0.0];
    gtr.unimodal_branch_likelihood = false;

    #[rustfmt::skip]
    let coefficients = array![
      [-0.04545042,  0.02261158,  0.02261158, 0.25],
      [ 0.04456328, -0.02228164, -0.02228164, 0.25],
    ];
    OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr))
  }

  /// Reconciliation against a positive optimizer output that is worse than
  /// zero on a multi-modal surface: the helper must return a positive mode
  /// (the better local max), not zero.
  ///
  /// Exercises the dispatch helper directly on the
  /// `positive_might_lose_to_zero` gate so that a refactor dropping the gate
  /// fails this test.
  ///
  /// The precondition asserts the three features of the Dinh-Matsen
  /// counterexample that the reconciliation relies on:
  ///
  /// - $\ell(0.2) > \ell(0)$: a positive point beats the boundary, so
  ///   the naive zero clamp would be wrong.
  /// - $\ell(0.2) > \ell(1.0)$: there is a local dip between the local
  ///   max and equilibrium, confirming the surface is multi-modal.
  /// - $\ell(5.0) > \ell(1.0)$: the function recovers from the dip
  ///   toward equilibrium, the recovery region the rustdoc on
  ///   `reconcile_zero_boundary` describes as "the global max at
  ///   equilibrium".
  #[test]
  fn test_dispatch_zero_boundary_reconcile_positive_candidate_finds_positive_mode() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let lh_zero = evaluate_mixed_log_lh_only(&contributions, 0.0);
    let lh_near_peak = evaluate_mixed_log_lh_only(&contributions, 0.2);
    let lh_at_dip = evaluate_mixed_log_lh_only(&contributions, 1.0);
    let lh_recovery = evaluate_mixed_log_lh_only(&contributions, 5.0);

    assert!(
      lh_near_peak > lh_zero,
      "precondition: log_lh(0.2)={lh_near_peak} > log_lh(0)={lh_zero}"
    );
    assert!(
      lh_near_peak > lh_at_dip,
      "precondition: local max log_lh(0.2)={lh_near_peak} > local dip log_lh(1.0)={lh_at_dip}"
    );
    assert!(
      lh_recovery > lh_at_dip,
      "precondition: recovery log_lh(5.0)={lh_recovery} > local dip log_lh(1.0)={lh_at_dip}"
    );

    let one_mutation = 0.01;
    // Simulate an inner solver returning a positive point in the local dip
    // at $t \approx 1$ where $\ell(\text{candidate}) < \ell(0)$. A naive
    // "force zero" would be wrong here because the local max at
    // $t \approx 0.2$ is better than both the dip and zero. The helper
    // must delegate to `grid_search_inner` and return the local max.
    let candidate = 1.0;
    let lh_candidate = evaluate_mixed_log_lh_only(&contributions, candidate);
    assert!(
      lh_candidate < lh_zero,
      "precondition: candidate log_lh({candidate})={lh_candidate} must be worse than log_lh(0)={lh_zero}"
    );

    // Grid extent covers $[0.1 \cdot \text{one\_mutation}, 0.5]$ because
    // the extent argument is smaller than the $0.5$ floor. That range
    // contains the local max at $t \approx 0.2$.
    let branch_length_extent = 0.2;
    let result =
      reconcile_zero_boundary(candidate, branch_length_extent, &contributions, 0, 0.0, one_mutation).unwrap();

    assert!(
      result > 0.0,
      "reconcile_zero_boundary must return a positive mode, not zero, got {result}"
    );
    let lh_result = evaluate_mixed_log_lh_only(&contributions, result);
    assert!(
      lh_result > lh_zero,
      "reconciled result must beat zero: log_lh(result)={lh_result} vs log_lh(0)={lh_zero}"
    );
    // The reconciled result lives in the grid range, and the local max
    // at $t \approx 0.2$ is the best positive grid candidate on that
    // range. Assert the result lies in a window around the local max.
    assert!(
      (0.1..0.4).contains(&result),
      "reconciled result should land near the local max at t ≈ 0.2, got {result}"
    );
  }

  /// Gate condition: `reconcile_zero_boundary` must pass a positive
  /// candidate through unchanged when any contribution has a degenerate
  /// site ($L_i(0) \leq 0$ or non-finite). The pass-through is enforced
  /// by `is_zero_better_than_grid_best`, which short-circuits on
  /// `all_sites_valid_at_zero` before calling the log-likelihood
  /// evaluator (which would produce $\ln(0) = -\infty$ or similar).
  /// Without this gate, a degenerate contribution would crash the
  /// evaluator or return nonsense.
  #[test]
  fn test_dispatch_zero_boundary_reconcile_degenerate_site_positive_candidate_passes_through() {
    // Degenerate JC69 site: $L_i(0) = 0$ so the log-likelihood at zero
    // is undefined. The contribution otherwise passes JC69's unimodal
    // classification.
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.0, 0.0, 0.0, 0.0]];
    let contribution = OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr));
    let contributions = [contribution];

    assert!(
      !contributions[0].all_sites_valid_at_zero(),
      "precondition: degenerate contribution must fail all_sites_valid_at_zero"
    );

    let candidate = 0.01;
    let result = reconcile_zero_boundary(candidate, 0.1, &contributions, 0, 0.0, 0.001).unwrap();
    assert!(
      result.to_bits() == candidate.to_bits(),
      "reconcile_zero_boundary must return positive candidate unchanged when a site is degenerate, got {result}"
    );
  }

  /// Gate condition: `reconcile_zero_boundary` must NOT pass an exact
  /// zero through when any contribution has a degenerate site. The
  /// run_optimize_mixed caller bumps a zero input branch length to
  /// `one_mutation` whenever `!all_sites_valid_at_zero`, precisely to
  /// keep the next marginal update in a well-defined evaluation domain.
  /// If the inner solver clamps back to zero, that result has the same
  /// problem and must be replaced with a positive grid candidate.
  ///
  /// Without this gate the helper would silently restore the invalid
  /// zero and the next marginal update would re-encounter $\ln(0)$ on
  /// the degenerate site.
  #[test]
  fn test_dispatch_zero_boundary_reconcile_degenerate_site_zero_candidate_routes_to_grid() {
    let gtr = jc69(JC69Params::default()).unwrap();
    // Two-site partition: site 0 is degenerate at $t = 0$ (all-zero
    // coefficients), site 1 is well-defined and prefers a positive
    // branch length.
    let coefficients = array![[0.0, 0.0, 0.0, 0.0], [0.5, 0.3, 0.1, 0.1]];
    let contribution = OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr));
    let contributions = [contribution];

    assert!(
      !contributions[0].all_sites_valid_at_zero(),
      "precondition: degenerate contribution must fail all_sites_valid_at_zero"
    );

    let candidate = 0.0;
    let one_mutation = 0.001;
    let result = reconcile_zero_boundary(candidate, 0.1, &contributions, 0, 0.0, one_mutation).unwrap();
    assert!(
      result > 0.0,
      "reconcile_zero_boundary must NOT preserve 0.0 when a site is degenerate at zero, got {result}"
    );
    assert!(result.is_finite(), "reconciled result must be finite, got {result}");
  }

  /// Gate condition: `reconcile_zero_boundary` must pass the candidate
  /// through unchanged when `indel_count > 0`. With indels present, the
  /// Poisson log-likelihood at $t = 0$ is $-\infty$ ($k \ln(\mu \cdot 0) - \mu \cdot 0 = k \cdot (-\infty) = -\infty$),
  /// so zero is never optimal and the boundary comparison must not fire.
  /// This is especially important because with a non-unimodal surface
  /// (Dinh-Matsen K80) the grid scan would otherwise run and find a
  /// positive mode, overriding a positive candidate that might be
  /// correct for the indel-weighted objective.
  #[test]
  fn test_dispatch_zero_boundary_reconcile_indel_count_positive_passes_through() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let candidate = 0.005;
    let indel_count = 1;
    let indel_rate = 44.4;
    let one_mutation = 0.01;

    let result =
      reconcile_zero_boundary(candidate, 0.2, &contributions, indel_count, indel_rate, one_mutation).unwrap();
    assert!(
      result.to_bits() == candidate.to_bits(),
      "reconcile_zero_boundary must return candidate unchanged when indel_count > 0, got {result}"
    );
  }

  /// End-to-end downstream check: after optimization on identical
  /// K80 sequences, `find_zero_optimal_internal_edges` must collect
  /// every internal edge as a collapse candidate. This closes the loop
  /// between "the post-dispatch reconciliation returns zero" and
  /// "topology cleanup actually removes the edge". The four methods
  /// that cannot evaluate at $t = 0$ (Brent, BrentSqrt, BrentLog,
  /// NewtonLog) place a strictly positive lower bound on the bracket
  /// (around $10^{-12}$) and rely on `reconcile_zero_boundary` to
  /// produce exact zero so that the `bl == 0.0` predicate in
  /// `find_zero_optimal_internal_edges` collects them.
  ///
  /// The 4-taxon rooted tree `((A,B)AB,(C,D)CD)root` has 2 internal
  /// edges (AB->root, CD->root) and 4 leaf edges (A,B,C,D). All 6
  /// should reach zero on identical sequences, but only the 2
  /// internal ones are eligible for collapse (leaves cannot be
  /// collapsed without removing their taxa).
  #[rustfmt::skip]
  #[rstest]
  #[case::brent(     BranchOptMethod::Brent)]
  #[case::brent_sqrt(BranchOptMethod::BrentSqrt)]
  #[case::brent_log( BranchOptMethod::BrentLog)]
  #[case::newton_log(BranchOptMethod::NewtonLog)]
  #[trace]
  fn test_dispatch_zero_boundary_topology_cleanup_collects_k80_internal_edges(
    #[case] method: BranchOptMethod,
  ) -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(IDENTICAL_TREE_NEWICK)?;
    let (_, sparse_partitions, mixed_partitions) = setup_identical_partitions(&graph, GtrModelName::K80)?;

    // Before optimization, every edge starts at 0.1 and none are zero.
    assert_eq!(
      0,
      find_zero_optimal_internal_edges(&graph, &sparse_partitions).len(),
      "precondition: no zero-length internal edges before optimization"
    );

    run_optimize_mixed(&graph, &mixed_partitions, method)?;

    let zero_edges = find_zero_optimal_internal_edges(&graph, &sparse_partitions);
    assert_eq!(
      2,
      zero_edges.len(),
      "{method:?}: expected 2 zero-optimal internal edges after optimization, got {}",
      zero_edges.len()
    );
    Ok(())
  }

  /// Pin the observed `newton_inner` behavior on the Dinh-Matsen K80
  /// surface: the $t$-space inner solver does NOT return exactly $0.0$
  /// from any of the 12 starting points spanning $t \in [0.05, 0.90]$
  /// on this multi-modal counterexample.
  ///
  /// # Context
  ///
  /// Both `newton_inner` and `newton_sqrt_inner` clamp the Newton
  /// step such that the iterate can hit exactly zero. On the same
  /// Dinh-Matsen K80 $\kappa = 3$ surface, `newton_sqrt_inner`
  /// actually does clamp to zero at $t_0 = 0.6$
  /// (see `test_dispatch_zero_boundary_newton_sqrt_inner_clamps_to_zero_on_dinh_matsen_k80`),
  /// and the exact-zero branch of `reconcile_zero_boundary` exists to
  /// catch that failure mode. This test pins the complementary
  /// observation: `newton_inner` in plain $t$-space is not affected
  /// because its chain rule produces a different local quadratic apex
  /// than the $s = \sqrt{t}$ reparameterization.
  ///
  /// # Theoretical analysis
  ///
  /// For Newton from `bl > 0` to clamp to zero, the local quadratic
  /// at `bl` (with derivative `d` and second derivative `d2 < 0`)
  /// must have its apex at `t* = bl - d/d2 <= 0`. In $t$-space on
  /// the Dinh-Matsen surface the apex stays positive at every
  /// starting point tested: Newton from $t_0 \in (0, 0.2)$ steps
  /// right toward the local max; Newton from $t_0 \in (0.2, 1)$
  /// steps left but the quadratic apex stays in the basin of the
  /// local max, not at zero. In $\sqrt{t}$-space the chain rule
  /// introduces extra factors of $2s$ and $4s^2$, and the resulting
  /// quadratic apex CAN hit zero at certain starting points.
  ///
  /// # What this test pins
  ///
  /// Across 12 starting points spanning the admissible interval on
  /// the Dinh-Matsen surface, `newton_inner` returns a strictly
  /// positive value. If a future refactor of `newton_inner` causes
  /// any starting point to start producing $0.0$ on this surface,
  /// the affected case fails and forces re-evaluation of whether
  /// `reconcile_zero_boundary` 's exact-zero branch needs to be
  /// extended or restructured.
  ///
  /// Starting points sweep the range that brackets the local maximum
  /// near $t = 0.2$, the local minimum near $t = 1.0$, and the
  /// recovery region toward equilibrium. They cover every concave
  /// segment of the multi-modal surface. Cases where the surface is
  /// non-concave at $t_0$ exercise `newton_inner`'s grid-search
  /// fallback rather than the clamping path; the positivity
  /// assertion holds for both since `grid_search_inner` returns a
  /// positive candidate when one beats zero.
  #[rustfmt::skip]
  #[rstest]
  #[case::t_0_05(0.05)]
  #[case::t_0_10(0.10)]
  #[case::t_0_15(0.15)]
  #[case::t_0_20(0.20)]
  #[case::t_0_25(0.25)]
  #[case::t_0_30(0.30)]
  #[case::t_0_40(0.40)]
  #[case::t_0_50(0.50)]
  #[case::t_0_60(0.60)]
  #[case::t_0_70(0.70)]
  #[case::t_0_80(0.80)]
  #[case::t_0_90(0.90)]
  #[trace]
  fn test_dispatch_zero_boundary_newton_inner_does_not_clamp_to_zero_on_dinh_matsen_k80(#[case] t0: f64) {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let one_mutation = 0.01;
    let metrics = evaluate_mixed(&contributions, t0);

    let result = newton_inner(t0, &metrics, &contributions, 0, 0.0, 0.0, one_mutation).unwrap();
    assert!(
      result > 0.0,
      "newton_inner from t0={t0} on the Dinh-Matsen K80 surface must return a positive value, got {result}. \
       If this assertion ever fails, `reconcile_zero_boundary` may need an exact-zero entry condition; \
       see the function rustdoc for the full justification."
    );
  }

  /// Positive reproduction: `newton_sqrt_inner` DOES clamp to exactly
  /// $0$ from $t_0 = 0.6$ on the Dinh-Matsen K80 counterexample,
  /// despite a better positive mode existing at $t \approx 0.2$.
  ///
  /// This is the structural step-clamping failure mode discussed in
  /// `reconcile_zero_boundary`'s rustdoc. The chain rule in $\sqrt{t}$
  /// space transforms the $t$-space derivatives to
  /// $d\ell/ds = 2s \cdot d\ell/dt$ and
  /// $d^2\ell/ds^2 = 4s^2 \cdot d^2\ell/dt^2 + 2 \cdot d\ell/dt$,
  /// which at $t = 0.6$ on this surface produces an $s$-space Hessian
  /// with a different sign than the $t$-space Hessian. The sqrt-space
  /// Newton step then exceeds the clamp bound `s`, and the clamp
  /// drives $s$ to $0$, giving $t = 0$. At that point the dispatch
  /// layer must recognize that the inner solver's zero return is NOT
  /// necessarily the global maximum and verify the positive domain.
  ///
  /// The exact-zero entry condition on `reconcile_zero_boundary`
  /// handles this reproduction: the full path is exercised by
  /// `test_dispatch_zero_boundary_reconcile_exact_zero_finds_positive_mode`.
  ///
  /// `newton_inner` does NOT exhibit this failure mode on the same
  /// surface (see `test_dispatch_zero_boundary_newton_inner_does_not_clamp_to_zero_on_dinh_matsen_k80`).
  /// The difference is entirely from the chain rule: the
  /// reparameterization changes the local quadratic apex.
  ///
  /// If a future refactor of `newton_sqrt_inner` changes the clamping
  /// behavior so that $t_0 = 0.6$ no longer clamps to $0$, this test
  /// fails and the rustdoc justification on `reconcile_zero_boundary`
  /// for the exact-zero gate must be revised.
  #[test]
  fn test_dispatch_zero_boundary_newton_sqrt_inner_clamps_to_zero_on_dinh_matsen_k80() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let one_mutation = 0.01;
    let t0 = 0.6;
    let metrics = evaluate_mixed(&contributions, t0);

    let result = newton_sqrt_inner(t0, &metrics, &contributions, 0, 0.0, 0.0, one_mutation).unwrap();
    assert!(
      result == 0.0,
      "reproduction: newton_sqrt_inner from t0={t0} on the Dinh-Matsen K80 surface must clamp to exactly 0 (else the reconcile zero-candidate gate is no longer necessary), got {result}"
    );
  }

  /// End-to-end for the exact-zero reconciliation path: given the
  /// inner solver returning $0$ on a multi-modal surface where a
  /// positive mode beats zero, `reconcile_zero_boundary` must reject
  /// the solver's output and return the positive mode via grid search.
  ///
  /// This is the downstream of
  /// `test_dispatch_zero_boundary_newton_sqrt_inner_clamps_to_zero_on_dinh_matsen_k80`:
  /// the sqrt solver returned $0$, the dispatch layer calls reconcile
  /// with `candidate = 0`, and the helper must detect the multi-modal
  /// context (non-unimodal model, all sites valid at zero, no indels)
  /// and fire grid search. Grid search finds the local max at
  /// $t \approx 0.2$ and returns it.
  #[test]
  fn test_dispatch_zero_boundary_reconcile_exact_zero_finds_positive_mode() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let lh_zero = evaluate_mixed_log_lh_only(&contributions, 0.0);
    let lh_near_peak = evaluate_mixed_log_lh_only(&contributions, 0.2);
    assert!(
      lh_near_peak > lh_zero,
      "precondition: log_lh(0.2)={lh_near_peak} > log_lh(0)={lh_zero}"
    );

    let one_mutation = 0.01;
    // Simulate newton_sqrt_inner returning exact zero by step clamping
    // from $t_0 = 0.6$ (reproduced by the companion pinning test).
    // The grid extent argument must still cover the local max at
    // $t \approx 0.2$; in practice `run_optimize_mixed` passes the
    // input branch length, which would be 0.6 in this scenario.
    let candidate = 0.0;
    let branch_length_extent = 0.6;
    let result =
      reconcile_zero_boundary(candidate, branch_length_extent, &contributions, 0, 0.0, one_mutation).unwrap();

    assert!(
      result > 0.0,
      "reconcile_zero_boundary must reject exact-zero and return a positive mode on a multi-modal surface, got {result}"
    );
    let lh_result = evaluate_mixed_log_lh_only(&contributions, result);
    assert!(
      lh_result > lh_zero,
      "reconciled result must beat zero: log_lh(result)={lh_result} vs log_lh(0)={lh_zero}"
    );
    assert!(
      (0.1..0.4).contains(&result),
      "reconciled result should land near the local max at t ≈ 0.2, got {result}"
    );
  }

  /// Exact-zero gate: on a unimodal model (JC69), an exact-zero
  /// candidate must pass through unchanged. The pre-dispatch shortcut
  /// would have fired if zero is the global max; otherwise the inner
  /// solver would have returned a positive mode. An exact-zero return
  /// reaching reconcile on an all-unimodal partition set is either
  /// correct (pass-through) or unreachable, so the grid scan is
  /// skipped to avoid wasted work.
  #[test]
  fn test_dispatch_zero_boundary_reconcile_exact_zero_unimodal_passes_through() {
    // JC69 with coefficients that assign weight to a non-zero eigenvalue,
    // so $\ell'(0) < 0$ and zero is the unique global max.
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.0, 1.0, 0.0, 0.0]];
    let contribution = OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr));
    let contributions = [contribution];

    assert!(
      contributions[0].has_unimodal_branch_likelihood(),
      "precondition: JC69 must be classified as unimodal"
    );

    let result = reconcile_zero_boundary(0.0, 0.1, &contributions, 0, 0.0, 0.01).unwrap();
    assert!(
      result == 0.0,
      "reconcile_zero_boundary must pass exact-zero through for unimodal models, got {result}"
    );
  }

  /// Exact-zero gate: with `indel_count > 0`, an exact-zero candidate
  /// must pass through unchanged regardless of model unimodality. With
  /// indels the Poisson log-likelihood at $t = 0$ is $-\infty$, so
  /// zero is never optimal; passing the candidate through is a no-op
  /// for correctness (the caller's initial branch length is clamped
  /// above zero for indel-bearing edges anyway).
  #[test]
  fn test_dispatch_zero_boundary_reconcile_exact_zero_indels_passes_through() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let result = reconcile_zero_boundary(0.0, 0.2, &contributions, 1, 44.4, 0.01).unwrap();
    assert!(
      result == 0.0,
      "reconcile_zero_boundary must pass exact-zero through when indels are present, got {result}"
    );
  }
}
