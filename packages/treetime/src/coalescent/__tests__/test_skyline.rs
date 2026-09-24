#[cfg(test)]
mod tests {
  use crate::coalescent::__tests__::helpers::{coalescent_node_times, constant_skyline, graph_with_dates};
  use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
  use crate::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use maplit::btreemap;
  use rstest::rstest;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_utils::assert_error;

  const SMALL_TREE_NWK: &str = "((leaf1:1.0,leaf2:1.0)internal1:1.0,leaf3:1.0)root:1.0;";

  fn small_tree_dates() -> DatesMap {
    scaled_small_tree_dates(1.0)
  }

  fn scaled_small_tree_dates(s: f64) -> DatesMap {
    let base = 2000.0;
    btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(base)),
      "internal1".to_owned() => Some(DateConstraint::exact(base + s * 5.0)),
      "leaf1".to_owned() => Some(DateConstraint::exact(base + s * 10.0)),
      "leaf2".to_owned() => Some(DateConstraint::exact(base + s * 10.0)),
      "leaf3".to_owned() => Some(DateConstraint::exact(base + s * 12.0)),
    }
  }

  #[test]
  fn test_optimize_skyline_returns_result() -> Result<(), Report> {
    let (graph, constraints) = graph_with_dates(SMALL_TREE_NWK, &small_tree_dates())?;
    let params = SkylineParams {
      n_points: 5,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params, &coalescent_node_times(&graph, &constraints))?;

    assert_eq!(6, result.segment_boundaries.len());
    assert_eq!(5, result.tc_schedule.values().len());
    pretty_assert_ulps_eq!(
      result.tc_schedule.breakpoints().view(),
      result.segment_boundaries.slice(ndarray::s![1..5]),
      max_ulps = 4
    );
    assert!(result.log_likelihood.value().is_finite());

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::negative(         -1.0,          "Skyline confidence must be finite and nonnegative, got -1")]
  #[case::not_a_number(     f64::NAN,      "Skyline confidence must be finite and nonnegative, got NaN")]
  #[case::positive_infinity(f64::INFINITY, "Skyline confidence must be finite and nonnegative, got inf")]
  #[trace]
  fn test_optimize_skyline_rejects_invalid_confidence(
    #[case] n_std: f64,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let (graph, constraints) = graph_with_dates(SMALL_TREE_NWK, &small_tree_dates())?;
    let params = SkylineParams {
      n_std,
      ..SkylineParams::default()
    };

    assert_error!(optimize_skyline(&graph, &params, &coalescent_node_times(&graph, &constraints)), expected);
    Ok(())
  }

  #[test]
  fn test_optimize_skyline_tc_values_positive() -> Result<(), Report> {
    let (graph, constraints) = graph_with_dates(SMALL_TREE_NWK, &small_tree_dates())?;
    let params = SkylineParams {
      n_points: 5,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params, &coalescent_node_times(&graph, &constraints))?;

    for &tc in result.tc_schedule.values() {
      assert!(
        tc > 0.0 && tc.is_finite(),
        "Tc segment value must be positive and finite, got {tc}"
      );
    }

    Ok(())
  }

  #[test]
  fn test_optimize_skyline_tc_distribution_evaluates() -> Result<(), Report> {
    let (graph, constraints) = graph_with_dates(SMALL_TREE_NWK, &small_tree_dates())?;
    let params = SkylineParams {
      n_points: 5,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params, &coalescent_node_times(&graph, &constraints))?;

    let t_min = result.segment_boundaries[0];
    let t_max = result.segment_boundaries[result.segment_boundaries.len() - 1];
    let t_mid = f64::midpoint(t_min, t_max);

    assert!(result.tc_distribution.eval(t_min)? > 0.0);
    assert!(result.tc_distribution.eval(t_mid)? > 0.0);
    assert!(result.tc_distribution.eval(t_max)? > 0.0);

    Ok(())
  }

  #[test]
  fn test_optimize_skyline_larger_tree() -> Result<(), Report> {
    const TREE_NWK: &str = "(((a:1,b:1)ab:1,(c:1,d:1)cd:1)abcd:1,((e:1,f:1)ef:1,(g:1,h:1)gh:1)efgh:1)root:1;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "abcd".to_owned() => Some(DateConstraint::exact(2002.0)),
      "efgh".to_owned() => Some(DateConstraint::exact(2002.0)),
      "ab".to_owned() => Some(DateConstraint::exact(2004.0)),
      "cd".to_owned() => Some(DateConstraint::exact(2004.0)),
      "ef".to_owned() => Some(DateConstraint::exact(2004.0)),
      "gh".to_owned() => Some(DateConstraint::exact(2004.0)),
      "a".to_owned() => Some(DateConstraint::exact(2006.0)),
      "b".to_owned() => Some(DateConstraint::exact(2007.0)),
      "c".to_owned() => Some(DateConstraint::exact(2008.0)),
      "d".to_owned() => Some(DateConstraint::exact(2009.0)),
      "e".to_owned() => Some(DateConstraint::exact(2010.0)),
      "f".to_owned() => Some(DateConstraint::exact(2011.0)),
      "g".to_owned() => Some(DateConstraint::exact(2012.0)),
      "h".to_owned() => Some(DateConstraint::exact(2013.0)),
    };

    let (graph, constraints) = graph_with_dates(TREE_NWK, &dates)?;
    let params = SkylineParams {
      n_points: 10,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params, &coalescent_node_times(&graph, &constraints))?;

    assert_eq!(10, result.tc_schedule.values().len());
    assert!(result.log_likelihood.value().is_finite());

    Ok(())
  }

  #[test]
  fn test_skyline_reported_likelihood_matches_model_evaluation() -> Result<(), Report> {
    let (graph, constraints) = graph_with_dates(SMALL_TREE_NWK, &small_tree_dates())?;
    let node_times = coalescent_node_times(&graph, &constraints);
    let params = SkylineParams {
      n_points: 4,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params, &node_times)?;
    let expected = compute_coalescent_total_lh(&graph, &result.tc_distribution, &node_times)?.value();

    pretty_assert_ulps_eq!(expected, result.log_likelihood.value(), max_ulps = 10);
    Ok(())
  }

  #[test]
  fn test_skyline_reported_likelihood_matches_per_edge_cost_for_polytomy() -> Result<(), Report> {
    const TREE_NWK: &str = "(a:1,b:1,c:1,d:1)root:1;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "a".to_owned() => Some(DateConstraint::exact(2010.0)),
      "b".to_owned() => Some(DateConstraint::exact(2010.0)),
      "c".to_owned() => Some(DateConstraint::exact(2010.0)),
      "d".to_owned() => Some(DateConstraint::exact(2010.0)),
    };
    let (graph, constraints) = graph_with_dates(TREE_NWK, &dates)?;
    let node_times = coalescent_node_times(&graph, &constraints);
    let params = SkylineParams {
      n_points: 2,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params, &node_times)?;
    let expected = compute_coalescent_total_lh(&graph, &result.tc_distribution, &node_times)?.value();

    pretty_assert_ulps_eq!(expected, result.log_likelihood.value(), max_ulps = 10);
    Ok(())
  }

  #[test]
  fn test_skyline_beats_or_matches_constant_tc() -> Result<(), Report> {
    let (graph, constraints) = graph_with_dates(SMALL_TREE_NWK, &small_tree_dates())?;
    let node_times = coalescent_node_times(&graph, &constraints);
    let params = SkylineParams {
      n_points: 4,
      stiffness: 1e-6,
      tolerance: 1e-12,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params, &node_times)?;

    let constant_tc = constant_skyline(&graph, &node_times)?;
    assert!(
      result.log_likelihood.value() >= constant_tc.log_likelihood.value() - 1e-10,
      "skyline LL {} should be >= constant-Tc LL {}",
      result.log_likelihood.value(),
      constant_tc.log_likelihood.value()
    );

    Ok(())
  }

  #[test]
  fn test_skyline_scale_invariant_trajectory() -> Result<(), Report> {
    let params = SkylineParams {
      n_points: 4,
      stiffness: 2.0,
      n_std: 2.0,
      tolerance: 1e-12,
      max_iter: 1000,
    };
    let s = 3.0;

    let (g1, c1) = graph_with_dates(SMALL_TREE_NWK, &scaled_small_tree_dates(1.0))?;
    let (gs, cs) = graph_with_dates(SMALL_TREE_NWK, &scaled_small_tree_dates(s))?;
    let r1 = optimize_skyline(&g1, &params, &coalescent_node_times(&g1, &c1))?;
    let rs = optimize_skyline(&gs, &params, &coalescent_node_times(&gs, &cs))?;

    for i in 0..params.n_points {
      let expected = s * r1.tc_schedule.values()[i];
      let rel = ((rs.tc_schedule.values()[i] - expected) / expected).abs();
      assert!(
        rel < 1e-10,
        "segment {i}: scaled Tc {} should be s×{} = {expected} (rel err {rel:.2e})",
        rs.tc_schedule.values()[i],
        r1.tc_schedule.values()[i]
      );
    }

    Ok(())
  }
}
