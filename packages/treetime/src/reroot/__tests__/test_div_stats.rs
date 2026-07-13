#[cfg(test)]
mod tests {
  use crate::payload::clock_set::ClockSet;
  use crate::reroot::div_stats::DivStats;
  use crate::reroot::traits::RootStats;
  use approx::assert_ulps_eq;

  // Oracle: divergence components of `ClockSet::leaf_contribution_to_parent`
  // (clock_set.rs:41-51), with `count` taking the role of `norm` and always
  // populated (no date gate).
  #[test]
  fn test_div_stats_leaf_contribution() {
    let stats = DivStats::leaf(None, 0.1, 1.0);
    assert_ulps_eq!(stats.count(), 1.0, max_ulps = 4);
    assert_ulps_eq!(stats.d_sum(), 0.1, max_ulps = 4);
    assert_ulps_eq!(stats.dsq_sum(), 0.01, max_ulps = 4);
  }

  #[test]
  fn test_div_stats_leaf_divides_by_variance() {
    let stats = DivStats::leaf(None, 0.2, 2.0);
    assert_ulps_eq!(stats.count(), 0.5, max_ulps = 4);
    assert_ulps_eq!(stats.d_sum(), 0.1, max_ulps = 4);
    assert_ulps_eq!(stats.dsq_sum(), 0.02, max_ulps = 4);
  }

  // Oracle: weighted variance of root-to-tip distances. For two tips at distances
  // d1, d2 with unit variance, score = (d1 - d2)^2 / 4 (analytical derivation);
  // also equals v0 base_regression chisq with fixed zero slope and no dates
  // (treeregression.py:45).
  #[test]
  fn test_div_stats_score_two_tips_analytical() {
    let d1 = 0.1;
    let d2 = 0.3;
    let stats = DivStats::new(2.0, d1 + d2, d1.powi(2) + d2.powi(2));
    let expected = (d1 - d2).powi(2) / 4.0;
    assert_ulps_eq!(stats.score(), expected, max_ulps = 4);
  }

  #[test]
  fn test_div_stats_score_zero_when_equidistant() {
    // Two tips at identical distances => zero variance => zero score.
    let stats = DivStats::new(2.0, 0.4, 2.0 * 0.2_f64.powi(2));
    assert_ulps_eq!(stats.score(), 0.0, max_ulps = 4);
  }

  // Oracle: `ClockSet::propagate_averages` (clock_set.rs:53-85). DivStats::propagate
  // must equal it on the divergence fields when all time terms are zero. A ClockSet
  // built from a tip with date 0 has t_sum = tsq_sum = dt_sum = 0 but norm = 1/var,
  // matching DivStats's `count`.
  #[test]
  fn test_div_stats_propagate_matches_clockset_with_zero_time() {
    let div = DivStats::leaf(None, 0.15, 1.5);
    let clock = ClockSet::leaf_contribution_to_parent(Some(0.0), 0.15, 1.5);

    let branch_length = 0.4;
    let variance = 0.7;
    let div_p = div.propagate(branch_length, variance);
    let clock_p = clock.propagate_averages(branch_length, variance);

    assert_ulps_eq!(div_p.count(), clock_p.norm(), max_ulps = 4);
    assert_ulps_eq!(div_p.d_sum(), clock_p.d_sum(), max_ulps = 4);
    assert_ulps_eq!(div_p.dsq_sum(), clock_p.dsq_sum(), max_ulps = 4);
    // The time terms on the ClockSet side stay identically zero.
    assert_ulps_eq!(clock_p.t_sum(), 0.0, max_ulps = 4);
    assert_ulps_eq!(clock_p.dt_sum(), 0.0, max_ulps = 4);
  }

  // With zero branch variance (the v0 internal-branch case), propagate reduces to
  // plain distance accumulation: count unchanged, distances shifted by bl.
  #[test]
  fn test_div_stats_propagate_zero_variance_accumulates() {
    let stats = DivStats::new(2.0, 0.4, 0.10);
    let bl = 0.3;
    let propagated = stats.propagate(bl, 0.0);
    assert_ulps_eq!(propagated.count(), 2.0, max_ulps = 4);
    assert_ulps_eq!(propagated.d_sum(), 0.4 + bl * 2.0, max_ulps = 4);
    assert_ulps_eq!(
      propagated.dsq_sum(),
      0.10 + 2.0 * bl * 0.4 + bl.powi(2) * 2.0,
      max_ulps = 4
    );
  }

  #[test]
  fn test_div_stats_add_is_elementwise() {
    let a = DivStats::new(1.0, 0.2, 0.04);
    let b = DivStats::new(2.0, 0.5, 0.13);
    let sum = a + b;
    assert_ulps_eq!(sum.count(), 3.0, max_ulps = 4);
    assert_ulps_eq!(sum.d_sum(), 0.7, max_ulps = 4);
    assert_ulps_eq!(sum.dsq_sum(), 0.17, max_ulps = 4);
  }

  #[test]
  fn test_div_stats_sub_inverts_add() {
    let a = DivStats::new(3.0, 0.7, 0.17);
    let b = DivStats::new(2.0, 0.5, 0.13);
    let diff = a + b - b;
    assert_ulps_eq!(diff.count(), a.count(), max_ulps = 4);
    assert_ulps_eq!(diff.d_sum(), a.d_sum(), max_ulps = 4);
    assert_ulps_eq!(diff.dsq_sum(), a.dsq_sum(), max_ulps = 4);
  }
}
