#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::commands::optimize::args::BranchOptMethod;
  use crate::commands::optimize::method_newton::newton_inner;
  use crate::commands::optimize::optimize_dense;
  use crate::commands::optimize::optimize_unified::{
    OptimizationContribution, evaluate_mixed, evaluate_mixed_log_lh_only, initial_guess_mixed, reconcile_zero_boundary,
    run_optimize_mixed,
  };
  use crate::commands::optimize::partition_ops::PartitionOptimizeVec;
  use crate::commands::optimize::run::collect_optimize_partitions;
  use crate::gtr::get_gtr::{JC69Params, K80Params, jc69, k80};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
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

  /// Set up dense + sparse partitions using K80 (non-unimodal) with identical
  /// leaf sequences. `initial_guess_mixed` is called with `overwrite_valid =
  /// false` so the positive branch lengths from the newick string survive
  /// into `run_optimize_mixed`, forcing each method to perform real
  /// optimization work.
  fn setup_k80_identical_partitions(graph: &GraphAncestral) -> Result<PartitionOptimizeVec, Report> {
    let aln = read_many_fasta_str(IDENTICAL_ALIGNMENT, &Alphabet::default())?;

    let dense_partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: k80(K80Params::default())?,
      alphabet: Alphabet::new(AlphabetName::Nuc)?,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    let sparse_partitions = vec![Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 1,
      gtr: k80(K80Params::default())?,
      alphabet: Alphabet::new(AlphabetName::Nuc)?,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(graph, &sparse_partitions, &aln)?;
    initialize_marginal(graph, &dense_partitions, &aln)?;
    update_marginal(graph, &sparse_partitions)?;

    let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);
    initial_guess_mixed(graph, &mixed_partitions, false)?;

    Ok(mixed_partitions)
  }

  /// All 6 methods must return exactly $t = 0$ on every edge for identical
  /// sequences under K80.
  ///
  /// K80 is not proven unimodal (Dinh & Matsen 2017, Corollary 3.1 does not
  /// apply), so the pre-dispatch `is_zero_branch_optimal` shortcut is
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
    let mixed_partitions = setup_k80_identical_partitions(&graph)?;

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
  /// This is the same scenario the original boundary check was designed to
  /// catch, but it now exercises the dispatch helper directly so that any
  /// future refactor that drops the gate gets caught.
  #[test]
  fn test_dispatch_zero_boundary_reconcile_positive_candidate_finds_positive_mode() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let lh_zero = evaluate_mixed_log_lh_only(&contributions, 0.0);
    let lh_positive_mode = evaluate_mixed_log_lh_only(&contributions, 0.2);
    assert!(
      lh_positive_mode > lh_zero,
      "precondition: K80 counterexample must satisfy log_lh(0.2)={lh_positive_mode} > log_lh(0)={lh_zero}"
    );

    let one_mutation = 0.01;
    // Simulate Brent or NewtonLog returning a tiny positive point that is
    // worse than zero (the bracket lower bound is strictly positive, so
    // these methods cannot evaluate exactly at zero).
    let candidate = 0.005;
    let branch_length_extent = 0.2;
    let result = reconcile_zero_boundary(candidate, branch_length_extent, &contributions, 0, 0.0, one_mutation);

    assert!(
      result > 0.0,
      "reconcile_zero_boundary must return a positive mode, not zero, got {result}"
    );
    let lh_result = evaluate_mixed_log_lh_only(&contributions, result);
    assert!(
      lh_result > lh_zero,
      "reconciled result must beat zero: log_lh(result)={lh_result} vs log_lh(0)={lh_zero}"
    );
  }

  /// Pin the observed `newton_inner` behavior on the Dinh-Matsen K80
  /// surface: the inner solver does NOT return exactly $0.0$ from any
  /// reasonable starting point on this multi-modal counterexample,
  /// despite the structural step-clamping pattern in `newton_inner`
  /// that could in principle return zero.
  ///
  /// # Why this test exists
  ///
  /// `reconcile_zero_boundary` only handles the `candidate > 0.0` case
  /// (NewtonLog and Brent variants returning a tiny positive value
  /// below the bracket floor). It does not handle an exact-zero return
  /// from the inner solver. The Newton clamp
  /// `new_bl = (bl - clamp(d/d2, -1.0, bl)).max(min)` can equal $0$
  /// when `d/d2 >= bl` with `d2 < 0`, and the first iteration on
  /// `bl = 0` always returns $0$ because the upper clamp bound is
  /// itself zero. Whether this scenario is reachable on a real
  /// non-unimodal surface where the global maximum is positive is the
  /// open question this test pins.
  ///
  /// # Theoretical analysis
  ///
  /// For Newton from `bl > 0` to clamp to zero, the local quadratic at
  /// `bl` (with derivative `d` and second derivative `d2 < 0`) must
  /// have its apex at `t* = bl - d/d2 <= 0`. That apex condition with
  /// `d2 < 0` means the local quadratic is monotone decreasing on
  /// `[0, bl]` with maximum at the boundary. For a competing positive
  /// mode `t** > bl` to satisfy $\ell(t^{**}) > \ell(0)$, the function
  /// must rise again past `bl`, which requires `d2` to change sign on
  /// $(bl, t^{**})$. Newton has no way to discover that other basin
  /// from a single Taylor expansion at `bl`; it converges within the
  /// local concave region containing `bl`.
  ///
  /// On the Dinh-Matsen K80 $\kappa = 3$ counterexample (Section 5,
  /// eq 5.1-5.2 of the cited paper), $\ell'(0) > 0$: the function is
  /// increasing at the boundary, so $0$ is not even a local maximum
  /// on this surface. The local max sits at $t \approx 0.2$. From
  /// any starting point in $(0, 0.2)$ Newton steps right toward the
  /// local max; from a starting point in $(0.2, 1)$ (between the
  /// local max and the local min near $t \approx 1$) Newton steps
  /// left but the quadratic apex stays in the basin of the local
  /// max, not at zero.
  ///
  /// # What this test pins
  ///
  /// Across 12 starting points spanning the admissible interval on
  /// the Dinh-Matsen surface, `newton_inner` returns a strictly
  /// positive value. If a future refactor of `newton_inner` changes
  /// the clamping bounds, the damping rule, or the inner-loop
  /// termination such that any starting point starts producing $0.0$
  /// on this surface, the affected case fails and forces
  /// re-evaluation of the helper's `candidate > 0.0` gate: either
  /// the new solver behavior is wrong (the surface has a real
  /// positive mode at $\approx 0.2$), or `reconcile_zero_boundary`
  /// needs to grow an exact-zero entry condition.
  ///
  /// Starting points sweep the range that brackets the local maximum
  /// near $t = 0.2$, the local minimum near $t = 1.0$, and the
  /// recovery region toward equilibrium. They cover every concave
  /// segment of the multi-modal surface. Cases where the surface is
  /// non-concave at $t_0$ exercise `newton_inner`'s grid-search
  /// fallback rather than the clamping path; the positivity
  /// assertion holds for both since `grid_search_inner` returns a
  /// positive candidate when one beats zero.
  ///
  /// # Limitations
  ///
  /// This test pins behavior only on one counterexample. Surfaces
  /// where $\ell'(0) < 0$ (zero is a local maximum but not the
  /// global maximum) could in principle trigger the structural
  /// clamping pattern. Constructing such a surface from a real GTR
  /// model requires a non-unimodal model on data where the boundary
  /// derivative is negative, which has not been exhibited in the
  /// codebase. If a real example is later identified, both this
  /// test and `reconcile_zero_boundary` need updating.
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

    let result = newton_inner(t0, &metrics, &contributions, 0, 0.0, 0.0, one_mutation);
    assert!(
      result > 0.0,
      "newton_inner from t0={t0} on the Dinh-Matsen K80 surface must return a positive value, got {result}. \
       If this assertion ever fails, `reconcile_zero_boundary` may need an exact-zero entry condition; \
       see the function rustdoc for the full justification."
    );
  }
}
