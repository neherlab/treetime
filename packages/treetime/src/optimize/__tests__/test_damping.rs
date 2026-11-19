#[cfg(test)]
mod tests {
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, compute_total_lh, setup_reconstruction, simple_alignment,
  };
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::iteration::apply_damping;
  use crate::optimize::loop_history::ConvergenceReason;
  use crate::optimize::params::OptimizeLoopSettings;
  use crate::optimize::params::{BranchOptMethod, TopologyOps};
  use crate::optimize::partition_rates::PartitionRates;
  use crate::optimize::run_loop::run_optimize_loop;
  use crate::progress::NoopProgress;
  use eyre::Report;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::key_map::EdgeMap;

  use treetime_io::nwk::nwk_read;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[test]
  fn test_apply_damping_zero_keeps_optimized_lengths() -> Result<(), Report> {
    let original = nwk_read(TREE_NEWICK.as_bytes())?.branch_lengths;
    let mut optimized: EdgeMap<Option<f64>> = original.iter().map(|(key, &bl)| (key, bl.map(|b| b * 2.0))).collect();
    let expected = optimized.clone();

    apply_damping(&mut optimized, &original, 0.0, 0);

    assert_eq!(expected, optimized);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::iter_0( 0, 0.750)]
  #[case::iter_1( 1, 0.5625)]
  #[case::iter_2( 2, 0.421875)]
  #[case::iter_4( 4, 0.2373046875)]
  #[case::iter_9( 9, 0.056313514709472656)]
  #[trace]
  fn test_apply_damping_weights_match_v0(#[case] iteration: usize, #[case] expected_old_weight: f64) -> Result<(), Report> {
    let old_bls = nwk_read(b"(A:1.0,B:1.0)root:0.0;".as_slice())?.branch_lengths;
    let mut bls: EdgeMap<Option<f64>> = old_bls.keys().map(|key| (key, Some(0.0))).collect();

    apply_damping(&mut bls, &old_bls, 0.75, iteration);

    let actual: Array1<f64> = bls.values().map(|bl| bl.expect("damped length")).collect();
    pretty_assert_abs_diff_eq!(array![expected_old_weight, expected_old_weight], actual, epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_apply_damping_blends_each_edge_with_its_own_old_length() -> Result<(), Report> {
    let old_bls = nwk_read(b"(A:0.1,B:0.2)root:0.0;".as_slice())?.branch_lengths;
    let mut bls: EdgeMap<Option<f64>> = old_bls.iter().map(|(key, &bl)| (key, bl.map(|b| b * 3.0))).collect();

    apply_damping(&mut bls, &old_bls, 0.75, 0);

    let damped_to_old_ratio: Array1<f64> = old_bls
      .keys()
      .map(|key| bls[key].expect("damped length") / old_bls[key].expect("old length"))
      .collect();
    pretty_assert_abs_diff_eq!(array![1.5, 1.5], damped_to_old_ratio, epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(BranchOptMethod::Newton)]
  #[case::brent( BranchOptMethod::Brent)]
  #[trace]
  fn test_damped_optimization_converges(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;
    let (reconstruction, initial_lh) = compute_total_lh(&graph, reconstruction, &branch_lengths)?;

    let max_iter = 10;
    let damping = 0.75;
    let dp = 0.1;

    let names_tt_2 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      &OptimizeLoopSettings {
        max_iter,
        dp,
        damping,
        opt_method: method,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
      branch_lengths,
      names_tt_2,
      &NoopProgress,
    )?;
    let reconstruction = result.partitions;

    assert!(
      matches!(
        result.stopped_at,
        Some((_, ConvergenceReason::Converged | ConvergenceReason::Oscillating))
      ),
      "Damped optimization did not converge within {max_iter} iterations: {:?}",
      result.stopped_at
    );

    let (_, final_lh) = compute_total_lh(&graph, reconstruction, &result.branch_lengths)?;
    assert!(
      final_lh > initial_lh,
      "Damped optimization did not improve the likelihood: {initial_lh:.6} -> {final_lh:.6}"
    );

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(BranchOptMethod::Newton)]
  #[case::brent( BranchOptMethod::Brent)]
  #[trace]
  fn test_damped_optimization_does_not_regress(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let (reconstruction, initial_lh) = compute_total_lh(&graph, reconstruction, &branch_lengths)?;

    let dp = 0.0;
    let names_tt_1 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      reconstruction,
      &OptimizeLoopSettings {
        max_iter: 10,
        dp,
        damping: 0.75,
        opt_method: method,
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
      branch_lengths,
      names_tt_1,
      &NoopProgress,
    )?;
    let reconstruction = result.partitions;

    let (_, final_lh) = compute_total_lh(&graph, reconstruction, &result.branch_lengths)?;
    assert!(
      final_lh >= initial_lh,
      "Damped optimization regressed: {initial_lh:.6} -> {final_lh:.6}"
    );
    Ok(())
  }
}
