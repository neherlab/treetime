#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_partitions, simple_alignment,
  };
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::method_brent::{brent_bracket, brent_log_inner, brent_sqrt_inner};
  use crate::optimize::method_newton::{chain_rule_log, chain_rule_sqrt};
  use crate::optimize::indel::{estimate_indel_rate, poisson_indel_log_lh};
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::likelihood::{OptimizationMetrics, evaluate_mixed, evaluate_mixed_log_lh_only, evaluate_with_indels_log_lh_only};
  use crate::optimize::method_newton::newton_tolerance_t;
  use crate::optimize::zero_boundary::min_branch_length_for_indels;
  use crate::partition::optimization_contribution::OptimizationContribution;
  use crate::partition::optimize_dense;
  use crate::partition::traits::PartitionOptimizeOps;
  use crate::partition::payload::ancestral::GraphAncestral;
  use crate::seq::indel::InDel;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use helpers::*;
  use ndarray::array;
  use parking_lot::RwLock;
  use rstest::rstest;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;

  /// At s=0, the first derivative is zero (chain rule factor 2s = 0) and
  /// the second derivative equals 2 * dl_dt (the only surviving term).
  #[test]
  fn test_optimize_method_chain_rule_at_zero() {
    let (ds, d2s) = chain_rule_sqrt(0.0, 100.0, -500.0);
    assert_abs_diff_eq!(ds, 0.0, epsilon = 1e-15);
    // d2l_ds2 = 4*0*(-500) + 2*100 = 200
    assert_abs_diff_eq!(d2s, 200.0, epsilon = 1e-15);
  }

  /// Chain rule transform: known analytical values.
  ///
  /// s=0.3 (t=0.09), dl_dt=10.0, d2l_dt2=-100.0:
  ///   dl_ds = 2 * 0.3 * 10.0 = 6.0
  ///   d2l_ds2 = 4 * 0.09 * (-100.0) + 2 * 10.0 = -36 + 20 = -16
  #[test]
  fn test_optimize_method_chain_rule_analytical() {
    let s = 0.3;
    let (ds, d2s) = chain_rule_sqrt(s, 10.0, -100.0);
    assert_abs_diff_eq!(ds, 6.0, epsilon = 1e-14);
    assert_abs_diff_eq!(d2s, -16.0, epsilon = 1e-14);
  }

  /// Chain rule log transform: known analytical values.
  ///
  /// t=0.09, dl_dt=10.0, d2l_dt2=-100.0:
  ///   dl_du = 0.09 * 10.0 = 0.9
  ///   d2l_du2 = 0.09^2 * (-100.0) + 0.09 * 10.0 = -0.81 + 0.9 = 0.09
  #[test]
  fn test_optimize_method_chain_rule_log_analytical() {
    let t = 0.09;
    let (du, d2u) = chain_rule_log(t, 10.0, -100.0);
    assert_abs_diff_eq!(du, 0.9, epsilon = 1e-14);
    assert_abs_diff_eq!(d2u, 0.09, epsilon = 1e-14);
  }

  /// At very small t, both log-space derivatives approach zero because
  /// the t factor suppresses them.
  #[test]
  fn test_optimize_method_chain_rule_log_small_t() {
    let t = 1e-10;
    let (du, d2u) = chain_rule_log(t, 1e6, -1e12);
    // dl_du = 1e-10 * 1e6 = 1e-4
    assert_abs_diff_eq!(du, 1e-4, epsilon = 1e-14);
    // d2l_du2 = (1e-10)^2 * (-1e12) + 1e-10 * 1e6 = -1e-8 + 1e-4 ≈ 1e-4
    assert_abs_diff_eq!(d2u, -1e-8 + 1e-4, epsilon = 1e-14);
  }

  /// u-space first derivative matches numerical central difference of the
  /// combined (substitution + indel) log-likelihood evaluated at t = exp(u).
  #[rustfmt::skip]
  #[rstest]
  #[case::small(   0.01,  2, 10.0)]
  #[case::medium(  0.1,   3, 20.0)]
  #[case::large(   0.5,   1,  5.0)]
  #[case::no_indel(0.2,   0,  0.0)]
  #[trace]
  fn test_optimize_method_chain_rule_log_numerical_first_derivative(
    #[case] t: f64,
    #[case] k: usize,
    #[case] mu: f64,
  ) {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.5, 0.3, 0.1, 0.1]];
    let contribution = OptimizationContribution::Dense(
      optimize_dense::PartitionContribution::new(coefficients, gtr),
    );
    let contributions = vec![contribution];

    let u = t.ln();
    let metrics = evaluate_mixed(&contributions, t);
    let indel = poisson_indel_log_lh(k, mu, t);
    let dl_dt = metrics.derivative + indel.derivative;
    let d2l_dt2 = metrics.second_derivative + indel.second_derivative;
    let (dl_du_analytical, _) = chain_rule_log(t, dl_dt, d2l_dt2);

    let h = u.abs() * 1e-5;
    let eval_u = |uv: f64| {
      let tv = uv.exp();
      evaluate_mixed_log_lh_only(&contributions, tv) + poisson_indel_log_lh(k, mu, tv).log_lh
    };
    let dl_du_numerical = (eval_u(u + h) - eval_u(u - h)) / (2.0 * h);

    // Central-difference first derivative: leading O(h^2) truncation at h = |u| * 1e-5
    // dominates round-off ~ eps_machine / h. 1e-4 is the tightest tolerance that holds
    // across the indel-heavy cases in this matrix.
    assert_abs_diff_eq!(dl_du_analytical, dl_du_numerical, epsilon = 1e-4);
  }

  /// u-space second derivative matches numerical central difference.
  #[rustfmt::skip]
  #[rstest]
  #[case::small(   0.01,  2, 10.0)]
  #[case::medium(  0.1,   3, 20.0)]
  #[case::large(   0.5,   1,  5.0)]
  #[case::no_indel(0.2,   0,  0.0)]
  #[trace]
  fn test_optimize_method_chain_rule_log_numerical_second_derivative(
    #[case] t: f64,
    #[case] k: usize,
    #[case] mu: f64,
  ) {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.5, 0.3, 0.1, 0.1]];
    let contribution = OptimizationContribution::Dense(
      optimize_dense::PartitionContribution::new(coefficients, gtr),
    );
    let contributions = vec![contribution];

    let u = t.ln();
    let metrics = evaluate_mixed(&contributions, t);
    let indel = poisson_indel_log_lh(k, mu, t);
    let dl_dt = metrics.derivative + indel.derivative;
    let d2l_dt2 = metrics.second_derivative + indel.second_derivative;
    let (_, d2l_du2_analytical) = chain_rule_log(t, dl_dt, d2l_dt2);

    let h = u.abs() * 1e-4;
    let eval_u = |uv: f64| {
      let tv = uv.exp();
      evaluate_mixed_log_lh_only(&contributions, tv) + poisson_indel_log_lh(k, mu, tv).log_lh
    };
    let d2l_du2_numerical = (eval_u(u + h) - 2.0 * eval_u(u) + eval_u(u - h)) / (h * h);

    // Central-difference second derivative: round-off ~ 4 * eps_machine * |f| / h^2.
    // At h = |u| * 1e-4 ~ 1e-5 with indel Hessian magnitude ~ 2e4 (k / t^2 at t=0.01,
    // k=2), the round-off bound is ~ 4 * 2.2e-16 * 2e4 / 1e-10 ~ 4e-2. 1e-2 is the
    // tightest tolerance that holds across the matrix; loosening further would mask
    // a real numerical defect.
    assert_abs_diff_eq!(d2l_du2_analytical, d2l_du2_numerical, epsilon = 1e-2);
  }

  /// s-space first derivative matches numerical central difference of the
  /// combined (substitution + indel) log-likelihood evaluated at t = s^2.
  #[rustfmt::skip]
  #[rstest]
  #[case::small(   0.01,  2, 10.0)]
  #[case::medium(  0.1,   3, 20.0)]
  #[case::large(   0.5,   1,  5.0)]
  #[case::no_indel(0.2,   0,  0.0)]
  #[trace]
  fn test_optimize_method_chain_rule_numerical_first_derivative(
    #[case] t: f64,
    #[case] k: usize,
    #[case] mu: f64,
  ) {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.5, 0.3, 0.1, 0.1]];
    let contribution = OptimizationContribution::Dense(
      optimize_dense::PartitionContribution::new(coefficients, gtr),
    );
    let contributions = vec![contribution];

    let s = t.sqrt();
    let metrics = evaluate_mixed(&contributions, t);
    let indel = poisson_indel_log_lh(k, mu, t);
    let dl_dt = metrics.derivative + indel.derivative;
    let d2l_dt2 = metrics.second_derivative + indel.second_derivative;
    let (dl_ds_analytical, _) = chain_rule_sqrt(s, dl_dt, d2l_dt2);

    let h = s * 1e-5;
    let eval_s = |sv: f64| {
      let tv = sv * sv;
      evaluate_mixed_log_lh_only(&contributions, tv) + poisson_indel_log_lh(k, mu, tv).log_lh
    };
    let dl_ds_numerical = (eval_s(s + h) - eval_s(s - h)) / (2.0 * h);

    // Central-difference first derivative: leading O(h^2) truncation at h = s * 1e-5
    // dominates round-off ~ eps_machine / h. 1e-4 is the tightest tolerance that holds
    // across the indel-heavy cases in this matrix.
    assert_abs_diff_eq!(dl_ds_analytical, dl_ds_numerical, epsilon = 1e-4);
  }

  /// s-space second derivative matches numerical central difference.
  #[rustfmt::skip]
  #[rstest]
  #[case::small(   0.01,  2, 10.0)]
  #[case::medium(  0.1,   3, 20.0)]
  #[case::large(   0.5,   1,  5.0)]
  #[case::no_indel(0.2,   0,  0.0)]
  #[trace]
  fn test_optimize_method_chain_rule_numerical_second_derivative(
    #[case] t: f64,
    #[case] k: usize,
    #[case] mu: f64,
  ) {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.5, 0.3, 0.1, 0.1]];
    let contribution = OptimizationContribution::Dense(
      optimize_dense::PartitionContribution::new(coefficients, gtr),
    );
    let contributions = vec![contribution];

    let s = t.sqrt();
    let metrics = evaluate_mixed(&contributions, t);
    let indel = poisson_indel_log_lh(k, mu, t);
    let dl_dt = metrics.derivative + indel.derivative;
    let d2l_dt2 = metrics.second_derivative + indel.second_derivative;
    let (_, d2l_ds2_analytical) = chain_rule_sqrt(s, dl_dt, d2l_dt2);

    let h = s * 1e-4;
    let eval_s = |sv: f64| {
      let tv = sv * sv;
      evaluate_mixed_log_lh_only(&contributions, tv) + poisson_indel_log_lh(k, mu, tv).log_lh
    };
    let d2l_ds2_numerical = (eval_s(s + h) - 2.0 * eval_s(s) + eval_s(s - h)) / (h * h);

    // Central-difference second derivative: round-off ~ 4 * eps_machine * |f| / h^2.
    // At h = s * 1e-4 ~ 1e-5 with indel Hessian magnitude ~ 2e4 (k / t^2 at t=0.01,
    // k=2), the round-off bound is ~ 4 * 2.2e-16 * 2e4 / 1e-10 ~ 4e-2. 1e-2 is the
    // tightest tolerance that holds across the matrix; loosening further would mask
    // a real numerical defect.
    assert_abs_diff_eq!(d2l_ds2_analytical, d2l_ds2_numerical, epsilon = 1e-2);
  }

  /// All methods produce finite non-negative branch lengths on a
  /// simple tree with no indels (well-conditioned substitution-only objective).
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
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let aln = simple_alignment()?;
    let (_, _, mixed_partitions) = setup_partitions(&graph, &aln)?;

    run_optimize_mixed(&graph, &mixed_partitions, method)?;

    for (i, edge_ref) in graph.get_edges().iter().enumerate() {
      let bl = edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(f64::NAN);
      assert!(bl.is_finite(), "Edge {i}: branch length is not finite ({bl})");
      assert!(bl >= 0.0, "Edge {i}: branch length is negative ({bl})");
    }

    Ok(())
  }

  /// C1: Local optimality. The combined log-likelihood at the reported
  /// optimum must exceed the log-likelihood at nearby points.
  ///
  /// Uses the indel rate captured before optimization (same rate the
  /// optimizer used) for consistent evaluation.
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
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (mixed_partitions, indel_rate) = setup_with_indels(&graph, 4)?;

    run_optimize_mixed(&graph, &mixed_partitions, method)?;

    let bl = first_edge_bl(&graph);
    assert!(bl > 0.0 && bl.is_finite(), "Optimized BL must be positive and finite, got {bl}");

    let lh_opt = eval_combined_first_edge(&graph, &mixed_partitions, indel_rate, bl)?;

    for &frac in &[0.001, 0.01, 0.1] {
      let delta = bl * frac;
      if bl - delta > 0.0 {
        let lh_below = eval_combined_first_edge(&graph, &mixed_partitions, indel_rate, bl - delta)?;
        assert!(
          lh_opt >= lh_below - 1e-10,
          "{method:?}: lh at t*={bl} ({lh_opt}) < lh at t*-{frac}*t ({lh_below})"
        );
      }
      let lh_above = eval_combined_first_edge(&graph, &mixed_partitions, indel_rate, bl + delta)?;
      assert!(
        lh_opt >= lh_above - 1e-10,
        "{method:?}: lh at t*={bl} ({lh_opt}) < lh at t*+{frac}*t ({lh_above})"
      );
    }

    Ok(())
  }

  /// C2: Stationarity. The implied Newton step at the optimum should be
  /// smaller than the Newton tolerance. Uses the optimizer's indel rate.
  #[rustfmt::skip]
  #[rstest]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[trace]
  fn test_optimize_method_stationarity(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (mixed_partitions, indel_rate) = setup_with_indels(&graph, 2)?;

    run_optimize_mixed(&graph, &mixed_partitions, method)?;

    let bl = first_edge_bl(&graph);
    let metrics = eval_metrics_first_edge(&graph, &mixed_partitions, indel_rate, bl)?;

    if metrics.second_derivative < 0.0 {
      let implied_step = (metrics.derivative / metrics.second_derivative).abs();
      let tol = newton_tolerance_t(bl);
      // Allow 10x tolerance: the optimizer stops when the step is below
      // tolerance, but the next step from the final position can be slightly
      // larger due to nonlinearity of the objective.
      assert!(
        implied_step < tol * 10.0,
        "{method:?}: implied Newton step ({implied_step}) exceeds 10x tolerance ({tol}), \
         dl={}, d2l={}, t*={bl}",
        metrics.derivative, metrics.second_derivative
      );
    }

    Ok(())
  }

  /// C3: Cross-method agreement. NewtonSqrt and Brent achieve similar
  /// combined log-likelihood values. Compares log-likelihood (not branch
  /// lengths) because the objective can be flat near the optimum.
  #[rustfmt::skip]
  #[rstest]
  #[case::k1(1)]
  #[case::k2(2)]
  #[case::k4(4)]
  #[trace]
  fn test_optimize_method_cross_method_lh_agreement(#[case] n_indels: usize) -> Result<(), Report> {
    let graph_brent: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (partitions_brent, rate_brent) = setup_with_indels(&graph_brent, n_indels)?;
    run_optimize_mixed(&graph_brent, &partitions_brent, BranchOptMethod::Brent)?;
    let bl_brent = first_edge_bl(&graph_brent);
    let lh_brent = eval_combined_first_edge(&graph_brent, &partitions_brent, rate_brent, bl_brent)?;

    let graph_sqrt: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (partitions_sqrt, rate_sqrt) = setup_with_indels(&graph_sqrt, n_indels)?;
    run_optimize_mixed(&graph_sqrt, &partitions_sqrt, BranchOptMethod::NewtonSqrt)?;
    let bl_sqrt = first_edge_bl(&graph_sqrt);
    let lh_sqrt = eval_combined_first_edge(&graph_sqrt, &partitions_sqrt, rate_sqrt, bl_sqrt)?;

    let lh_diff = (lh_brent - lh_sqrt).abs();
    assert!(
      lh_diff < 1e-3,
      "Brent lh ({lh_brent}) and NewtonSqrt lh ({lh_sqrt}) differ by {lh_diff}, \
       bl_brent={bl_brent}, bl_sqrt={bl_sqrt}"
    );

    Ok(())
  }

  /// C3b: Cross-method agreement. NewtonLog and Brent achieve similar
  /// combined log-likelihood values.
  #[rustfmt::skip]
  #[rstest]
  #[case::k1(1)]
  #[case::k2(2)]
  #[case::k4(4)]
  #[trace]
  fn test_optimize_method_cross_method_lh_agreement_newton_log(#[case] n_indels: usize) -> Result<(), Report> {
    let graph_brent: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (partitions_brent, rate_brent) = setup_with_indels(&graph_brent, n_indels)?;
    run_optimize_mixed(&graph_brent, &partitions_brent, BranchOptMethod::Brent)?;
    let bl_brent = first_edge_bl(&graph_brent);
    let lh_brent = eval_combined_first_edge(&graph_brent, &partitions_brent, rate_brent, bl_brent)?;

    let graph_log: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (partitions_log, rate_log) = setup_with_indels(&graph_log, n_indels)?;
    run_optimize_mixed(&graph_log, &partitions_log, BranchOptMethod::NewtonLog)?;
    let bl_log = first_edge_bl(&graph_log);
    let lh_log = eval_combined_first_edge(&graph_log, &partitions_log, rate_log, bl_log)?;

    let lh_diff = (lh_brent - lh_log).abs();
    assert!(
      lh_diff < 1e-3,
      "Brent lh ({lh_brent}) and NewtonLog lh ({lh_log}) differ by {lh_diff}, \
       bl_brent={bl_brent}, bl_log={bl_log}"
    );

    Ok(())
  }

  /// C3c: Full cross-method log-likelihood agreement.
  ///
  /// Each non-reference method's combined log-likelihood must agree with
  /// `BrentSqrt`'s within 1e-3 on the same input. Compares log-likelihood
  /// (not branch lengths) because the objective can be flat near the optimum.
  ///
  /// Parameterized one case per (method, n_indels) so a regression in any
  /// single method on any single indel count fails its own case rather than
  /// short-circuiting on the first failure across a manual loop.
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
    // Reference: run BrentSqrt (the v0-matching default) on a fresh graph
    // and capture its post-optimization log-likelihood at the first edge.
    let lh_ref = {
      let graph_ref: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
      let (partitions_ref, rate_ref) = setup_with_indels(&graph_ref, n_indels)?;
      run_optimize_mixed(&graph_ref, &partitions_ref, BranchOptMethod::BrentSqrt)?;
      let bl_ref = first_edge_bl(&graph_ref);
      eval_combined_first_edge(&graph_ref, &partitions_ref, rate_ref, bl_ref)?
    };

    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (partitions, rate) = setup_with_indels(&graph, n_indels)?;
    run_optimize_mixed(&graph, &partitions, method)?;
    let bl = first_edge_bl(&graph);
    let lh = eval_combined_first_edge(&graph, &partitions, rate, bl)?;

    let diff = (lh - lh_ref).abs();
    assert!(
      diff < 1e-3,
      "{method:?} (n_indels={n_indels}) lh ({lh}) differs from BrentSqrt lh ({lh_ref}) by {diff} > 1e-3, bl={bl}"
    );

    Ok(())
  }

  /// C5b: NewtonLog achieves equal or better log-likelihood than Newton
  /// in t-space on the Hessian-dominated case.
  #[test]
  fn test_optimize_method_newton_log_improves_over_newton() -> Result<(), Report> {
    let graph_newton: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (partitions_newton, rate_newton) = setup_with_indels(&graph_newton, 4)?;
    run_optimize_mixed(&graph_newton, &partitions_newton, BranchOptMethod::Newton)?;
    let bl_newton = first_edge_bl(&graph_newton);
    let lh_newton = eval_combined_first_edge(&graph_newton, &partitions_newton, rate_newton, bl_newton)?;

    let graph_log: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (partitions_log, rate_log) = setup_with_indels(&graph_log, 4)?;
    run_optimize_mixed(&graph_log, &partitions_log, BranchOptMethod::NewtonLog)?;
    let bl_log = first_edge_bl(&graph_log);
    let lh_log = eval_combined_first_edge(&graph_log, &partitions_log, rate_log, bl_log)?;

    assert!(bl_newton > 0.0 && bl_newton.is_finite());
    assert!(bl_log > 0.0 && bl_log.is_finite());

    assert!(
      lh_log >= lh_newton - 1e-10,
      "NewtonLog lh ({lh_log}) should be >= Newton lh ({lh_newton}), \
       bl_log={bl_log}, bl_newton={bl_newton}"
    );

    Ok(())
  }

  /// C5: NewtonSqrt achieves equal or better log-likelihood than Newton
  /// in t-space on the Hessian-dominated case.
  #[test]
  fn test_optimize_method_newton_sqrt_improves_over_newton() -> Result<(), Report> {
    let graph_newton: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (partitions_newton, rate_newton) = setup_with_indels(&graph_newton, 4)?;
    run_optimize_mixed(&graph_newton, &partitions_newton, BranchOptMethod::Newton)?;
    let bl_newton = first_edge_bl(&graph_newton);
    let lh_newton = eval_combined_first_edge(&graph_newton, &partitions_newton, rate_newton, bl_newton)?;

    let graph_sqrt: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (partitions_sqrt, rate_sqrt) = setup_with_indels(&graph_sqrt, 4)?;
    run_optimize_mixed(&graph_sqrt, &partitions_sqrt, BranchOptMethod::NewtonSqrt)?;
    let bl_sqrt = first_edge_bl(&graph_sqrt);
    let lh_sqrt = eval_combined_first_edge(&graph_sqrt, &partitions_sqrt, rate_sqrt, bl_sqrt)?;

    assert!(bl_newton > 0.0 && bl_newton.is_finite());
    assert!(bl_sqrt > 0.0 && bl_sqrt.is_finite());

    assert!(
      lh_sqrt >= lh_newton - 1e-10,
      "NewtonSqrt lh ({lh_sqrt}) should be >= Newton lh ({lh_newton}), \
       bl_sqrt={bl_sqrt}, bl_newton={bl_newton}"
    );

    Ok(())
  }

  /// C5c: Cross-conditioning ordering for Newton variants on indel-bearing edges.
  ///
  /// Better-conditioned parameterizations find equal or better optima:
  /// `lh_newton_log >= lh_newton_sqrt >= lh_newton` within tolerance.
  ///
  /// This documents the known Newton-t limitation: it is a correct implementation
  /// of a limited algorithm, not a bug. If Newton-t produces better log-likelihood
  /// than Newton-log, the implementation is wrong.
  #[rustfmt::skip]
  #[rstest]
  #[case::k1(1)]
  #[case::k2(2)]
  #[case::k4(4)]
  #[trace]
  fn test_optimize_method_newton_cross_conditioning_ordering(#[case] n_indels: usize) -> Result<(), Report> {
    let graph_newton: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (partitions_newton, rate_newton) = setup_with_indels(&graph_newton, n_indels)?;
    run_optimize_mixed(&graph_newton, &partitions_newton, BranchOptMethod::Newton)?;
    let bl_newton = first_edge_bl(&graph_newton);
    let lh_newton = eval_combined_first_edge(&graph_newton, &partitions_newton, rate_newton, bl_newton)?;

    let graph_sqrt: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (partitions_sqrt, rate_sqrt) = setup_with_indels(&graph_sqrt, n_indels)?;
    run_optimize_mixed(&graph_sqrt, &partitions_sqrt, BranchOptMethod::NewtonSqrt)?;
    let bl_sqrt = first_edge_bl(&graph_sqrt);
    let lh_sqrt = eval_combined_first_edge(&graph_sqrt, &partitions_sqrt, rate_sqrt, bl_sqrt)?;

    let graph_log: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (partitions_log, rate_log) = setup_with_indels(&graph_log, n_indels)?;
    run_optimize_mixed(&graph_log, &partitions_log, BranchOptMethod::NewtonLog)?;
    let bl_log = first_edge_bl(&graph_log);
    let lh_log = eval_combined_first_edge(&graph_log, &partitions_log, rate_log, bl_log)?;

    // Verify ordering: lh_newton_log >= lh_newton_sqrt >= lh_newton
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

  /// All Brent variants produce positive, finite branch lengths with indels.
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
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (mixed_partitions, _) = setup_with_indels(&graph, n_indels)?;

    run_optimize_mixed(&graph, &mixed_partitions, method)?;

    let bl = first_edge_bl(&graph);
    assert!(bl > 0.0, "{method:?} BL with {n_indels} indels must be positive, got {bl}");
    assert!(bl.is_finite(), "{method:?} BL must be finite, got {bl}");
    Ok(())
  }

  /// NewtonLog produces positive, finite branch lengths with indels present.
  #[rustfmt::skip]
  #[rstest]
  #[case::k1(1)]
  #[case::k2(2)]
  #[case::k4(4)]
  #[trace]
  fn test_optimize_method_newton_log_positive_with_indels(#[case] n_indels: usize) -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (mixed_partitions, _) = setup_with_indels(&graph, n_indels)?;

    run_optimize_mixed(&graph, &mixed_partitions, BranchOptMethod::NewtonLog)?;

    let bl = first_edge_bl(&graph);
    assert!(bl > 0.0, "NewtonLog BL with {n_indels} indels must be positive, got {bl}");
    assert!(bl.is_finite(), "NewtonLog BL must be finite, got {bl}");
    Ok(())
  }

  /// NewtonSqrt produces positive, finite branch lengths with indels present.
  #[rustfmt::skip]
  #[rstest]
  #[case::k1(1)]
  #[case::k2(2)]
  #[case::k4(4)]
  #[trace]
  fn test_optimize_method_newton_sqrt_positive_with_indels(#[case] n_indels: usize) -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (mixed_partitions, _) = setup_with_indels(&graph, n_indels)?;

    run_optimize_mixed(&graph, &mixed_partitions, BranchOptMethod::NewtonSqrt)?;

    let bl = first_edge_bl(&graph);
    assert!(bl > 0.0, "NewtonSqrt BL with {n_indels} indels must be positive, got {bl}");
    assert!(bl.is_finite(), "NewtonSqrt BL must be finite, got {bl}");
    Ok(())
  }

  /// Newton (t-space) produces positive, finite branch lengths with indels present.
  #[rustfmt::skip]
  #[rstest]
  #[case::k1(1)]
  #[case::k2(2)]
  #[case::k4(4)]
  #[trace]
  fn test_optimize_method_newton_positive_with_indels(#[case] n_indels: usize) -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (mixed_partitions, _) = setup_with_indels(&graph, n_indels)?;

    run_optimize_mixed(&graph, &mixed_partitions, BranchOptMethod::Newton)?;

    let bl = first_edge_bl(&graph);
    assert!(bl > 0.0, "Newton BL with {n_indels} indels must be positive, got {bl}");
    assert!(bl.is_finite(), "Newton BL must be finite, got {bl}");
    Ok(())
  }

  /// All three Brent parameterizations achieve similar log-likelihood.
  #[test]
  fn test_optimize_method_brent_cross_parameterization_lh_agreement() -> Result<(), Report> {
    let n_indels = 3;

    let graph_t: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (parts_t, rate_t) = setup_with_indels(&graph_t, n_indels)?;
    run_optimize_mixed(&graph_t, &parts_t, BranchOptMethod::Brent)?;
    let bl_t = first_edge_bl(&graph_t);
    let lh_t = eval_combined_first_edge(&graph_t, &parts_t, rate_t, bl_t)?;

    let graph_sqrt: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (parts_sqrt, rate_sqrt) = setup_with_indels(&graph_sqrt, n_indels)?;
    run_optimize_mixed(&graph_sqrt, &parts_sqrt, BranchOptMethod::BrentSqrt)?;
    let bl_sqrt = first_edge_bl(&graph_sqrt);
    let lh_sqrt = eval_combined_first_edge(&graph_sqrt, &parts_sqrt, rate_sqrt, bl_sqrt)?;

    let graph_log: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (parts_log, rate_log) = setup_with_indels(&graph_log, n_indels)?;
    run_optimize_mixed(&graph_log, &parts_log, BranchOptMethod::BrentLog)?;
    let bl_log = first_edge_bl(&graph_log);
    let lh_log = eval_combined_first_edge(&graph_log, &parts_log, rate_log, bl_log)?;

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

  /// brent_sqrt_inner: the transform s=sqrt(t), evaluate at s^2, square
  /// result back must produce a branch length in t-space that is a local
  /// optimum of the original objective.
  #[test]
  fn test_optimize_method_brent_sqrt_transform_round_trip() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let aln = simple_alignment()?;
    let (_, _, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let edge_key = graph.get_edges()[0].read_arc().key();
    let contributions: Vec<OptimizationContribution> = mixed_partitions
      .iter()
      .map(|p| p.read_arc().create_edge_contribution(edge_key))
      .collect::<Result<_, _>>()?;

    let total_length: usize = mixed_partitions.iter().map(|p| p.read_arc().sequence_length()).sum();
    let one_mutation = 1.0 / total_length as f64;
    let branch_length = 0.01;

    let result = brent_sqrt_inner(branch_length, &contributions, 0, 0.0, 0.0, one_mutation).unwrap();
    assert!(
      result >= 0.0,
      "brent_sqrt_inner result must be non-negative, got {result}"
    );
    assert!(
      result.is_finite(),
      "brent_sqrt_inner result must be finite, got {result}"
    );

    // Verify local optimality in t-space
    let lh_opt = evaluate_with_indels_log_lh_only(&contributions, 0, 0.0, result);
    if result > 1e-10 {
      let lh_below = evaluate_with_indels_log_lh_only(&contributions, 0, 0.0, result * 0.99);
      let lh_above = evaluate_with_indels_log_lh_only(&contributions, 0, 0.0, result * 1.01);
      assert!(
        lh_opt >= lh_below - 1e-10,
        "sqrt: lh at opt ({lh_opt}) < lh below ({lh_below})"
      );
      assert!(
        lh_opt >= lh_above - 1e-10,
        "sqrt: lh at opt ({lh_opt}) < lh above ({lh_above})"
      );
    }

    Ok(())
  }

  /// brent_log_inner: the transform u=ln(t), evaluate at exp(u), exponentiate
  /// result back must produce a branch length in t-space that is a local
  /// optimum of the original objective.
  #[test]
  fn test_optimize_method_brent_log_transform_round_trip() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let aln = simple_alignment()?;
    let (_, _, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let edge_key = graph.get_edges()[0].read_arc().key();
    let contributions: Vec<OptimizationContribution> = mixed_partitions
      .iter()
      .map(|p| p.read_arc().create_edge_contribution(edge_key))
      .collect::<Result<_, _>>()?;

    let total_length: usize = mixed_partitions.iter().map(|p| p.read_arc().sequence_length()).sum();
    let one_mutation = 1.0 / total_length as f64;
    let branch_length = 0.01;

    let result = brent_log_inner(branch_length, &contributions, 0, 0.0, 0.0, one_mutation).unwrap();
    assert!(result > 0.0, "brent_log_inner result must be positive, got {result}");
    assert!(
      result.is_finite(),
      "brent_log_inner result must be finite, got {result}"
    );

    // Verify local optimality in t-space
    let lh_opt = evaluate_with_indels_log_lh_only(&contributions, 0, 0.0, result);
    let lh_below = evaluate_with_indels_log_lh_only(&contributions, 0, 0.0, result * 0.99);
    let lh_above = evaluate_with_indels_log_lh_only(&contributions, 0, 0.0, result * 1.01);
    assert!(
      lh_opt >= lh_below - 1e-10,
      "log: lh at opt ({lh_opt}) < lh below ({lh_below})"
    );
    assert!(
      lh_opt >= lh_above - 1e-10,
      "log: lh at opt ({lh_opt}) < lh above ({lh_above})"
    );

    Ok(())
  }

  /// C4: Brent bracket validity for all three parameterizations.
  ///
  /// The bracket used by production code is computed from the *input* branch
  /// length (before optimization), not the optimized result. This test captures
  /// the input BL first, computes the bracket from it, runs optimization, then
  /// verifies the optimum beats both bracket endpoints.
  #[rustfmt::skip]
  #[rstest]
  #[case::brent(     BranchOptMethod::Brent)]
  #[case::brent_sqrt(BranchOptMethod::BrentSqrt)]
  #[case::brent_log( BranchOptMethod::BrentLog)]
  #[trace]
  fn test_optimize_method_brent_bracket_validity(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (mixed_partitions, indel_rate) = setup_with_indels(&graph, 4)?;

    // Capture input branch length and compute bracket BEFORE optimization
    // (production code computes the bracket from the input BL). Calls the
    // production helpers brent_bracket and min_branch_length_for_indels
    // directly to avoid drift if either formula changes.
    let input_bl = graph.get_edges()[0].read_arc().payload().read_arc().branch_length().unwrap_or(0.0);
    let total_length: usize = mixed_partitions.iter().map(|p| p.read_arc().sequence_length()).sum();
    let one_mutation = 1.0 / total_length as f64;
    let min_bl = min_branch_length_for_indels(4, one_mutation);
    let (lower, upper) = brent_bracket(input_bl, min_bl, one_mutation);

    run_optimize_mixed(&graph, &mixed_partitions, method)?;

    let bl = first_edge_bl(&graph);
    let lh_opt = eval_combined_first_edge(&graph, &mixed_partitions, indel_rate, bl)?;

    let lh_lower = eval_combined_first_edge(&graph, &mixed_partitions, indel_rate, lower)?;
    let lh_upper = eval_combined_first_edge(&graph, &mixed_partitions, indel_rate, upper)?;

    assert!(
      lh_opt >= lh_lower - 1e-10,
      "{method:?} optimum lh ({lh_opt}) < lower bracket lh ({lh_lower})"
    );
    assert!(
      lh_opt >= lh_upper - 1e-10,
      "{method:?} optimum lh ({lh_opt}) < upper bracket lh ({lh_upper})"
    );

    Ok(())
  }

  mod generators {
    use proptest::prelude::*;
    /// Strategy for `s = sqrt(t)` values used by `chain_rule_sqrt`.
    /// Bounded above well below the catastrophic-cancellation regime of the
    /// step-clamping invariants (production never sees `s` near these caps).
    pub fn gen_s() -> impl Strategy<Value = f64> {
      1e-6_f64..1e3_f64
    }
    /// Strategy for `t > 0` values used by `chain_rule_log`.
    pub fn gen_t() -> impl Strategy<Value = f64> {
      1e-10_f64..1e3_f64
    }
    /// Strategy for first-derivative magnitudes typical of substitution +
    /// indel objectives at the per-edge scale.
    pub fn gen_dl_dt() -> impl Strategy<Value = f64> {
      -1e6_f64..1e6_f64
    }
    /// Strategy for second-derivative magnitudes (Hessians can swing widely
    /// on indel-heavy short branches).
    pub fn gen_d2l_dt2() -> impl Strategy<Value = f64> {
      -1e8_f64..1e8_f64
    }
    /// Strategy for a non-zero scalar multiplier in the linearity invariant.
    pub fn gen_scalar() -> impl Strategy<Value = f64> {
      prop_oneof![-1e3_f64..-1e-3_f64, 1e-3_f64..1e3_f64]
    }
  }

  use proptest::prelude::*;

  proptest! {
    /// `chain_rule_sqrt(s, dl_dt, d2l_dt2)` must equal the closed-form
    /// (2*s*dl_dt, 4*s^2*d2l_dt2 + 2*dl_dt). Pinning the formula identity
    /// catches any algebraic regression in the chain rule.
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

    /// `chain_rule_sqrt` is linear in `(dl_dt, d2l_dt2)`: scaling both inputs
    /// by `k` scales both outputs by `k`. This is the homogeneity invariant
    /// of a linear transform.
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

    /// `chain_rule_log(t, dl_dt, d2l_dt2)` must equal the closed-form
    /// (t*dl_dt, t^2*d2l_dt2 + t*dl_dt).
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

    /// `chain_rule_log` is linear in `(dl_dt, d2l_dt2)`.
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

  mod helpers {
    use super::*;

    /// Set up dense+sparse partitions from simple_alignment() (which has
    /// real substitutions between taxa) and inject indels on the first edge.
    /// Returns (mixed_partitions, indel_rate_before_optimization).
    ///
    /// The indel rate is captured BEFORE optimization because
    /// `run_optimize_mixed` computes the rate once from initial branch
    /// lengths and uses it as a constant throughout. Test evaluations must
    /// use the same rate for consistent verification.
    pub(super) fn setup_with_indels(
      graph: &GraphAncestral,
      n_indels: usize,
    ) -> Result<(crate::partition::traits::PartitionOptimizeVec, f64), Report> {
      let aln = simple_alignment()?;
      let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(graph, &aln)?;

      let first_edge_key = graph.get_edges()[0].read_arc().key();
      let indels: Vec<InDel> = (0..n_indels)
        .map(|i| InDel::del((i * 3, i * 3 + 3), Seq::try_from_str("ACG").unwrap()))
        .collect();

      for p in &dense_partitions {
        p.write_arc().data.edges.get_mut(&first_edge_key).unwrap().indels = indels.clone();
      }
      for p in &sparse_partitions {
        p.write_arc().edges.get_mut(&first_edge_key).unwrap().indels = indels.clone();
      }

      // Capture indel rate at the same point run_optimize_mixed will
      let indel_rate = estimate_indel_rate(graph, &mixed_partitions);

      Ok((mixed_partitions, indel_rate))
    }

    /// Evaluate combined (substitution + indel) log-likelihood at branch
    /// length t for the first edge, using a fixed indel rate (the rate
    /// the optimizer used, not the post-optimization rate).
    pub(super) fn eval_combined_first_edge(
      graph: &GraphAncestral,
      partitions: &[Arc<RwLock<dyn PartitionOptimizeOps>>],
      indel_rate: f64,
      t: f64,
    ) -> Result<f64, Report> {
      let edge_key = graph.get_edges()[0].read_arc().key();
      let contributions: Vec<OptimizationContribution> = partitions
        .iter()
        .map(|p| p.read_arc().create_edge_contribution(edge_key))
        .collect::<Result<_, _>>()?;

      let indel_count: usize = partitions.iter().map(|p| p.read_arc().edge_indel_count(edge_key)).sum();

      let sub_lh = evaluate_mixed_log_lh_only(&contributions, t);
      let indel_lh = poisson_indel_log_lh(indel_count, indel_rate, t).log_lh;
      Ok(sub_lh + indel_lh)
    }

    /// Evaluate combined metrics at branch length t for the first edge,
    /// using a fixed indel rate.
    pub(super) fn eval_metrics_first_edge(
      graph: &GraphAncestral,
      partitions: &[Arc<RwLock<dyn PartitionOptimizeOps>>],
      indel_rate: f64,
      t: f64,
    ) -> Result<OptimizationMetrics, Report> {
      let edge_key = graph.get_edges()[0].read_arc().key();
      let contributions: Vec<OptimizationContribution> = partitions
        .iter()
        .map(|p| p.read_arc().create_edge_contribution(edge_key))
        .collect::<Result<_, _>>()?;

      let indel_count: usize = partitions.iter().map(|p| p.read_arc().edge_indel_count(edge_key)).sum();

      let mut metrics = evaluate_mixed(&contributions, t);
      metrics.add(&poisson_indel_log_lh(indel_count, indel_rate, t));
      Ok(metrics)
    }

    /// Branch length on the first edge after optimization, panicking if it
    /// is missing or NaN. Captures the 5-line read chain that recurs across
    /// tests in this file.
    pub(super) fn first_edge_bl(graph: &GraphAncestral) -> f64 {
      graph.get_edges()[0]
        .read_arc()
        .payload()
        .read_arc()
        .branch_length()
        .unwrap()
    }
  }
}
