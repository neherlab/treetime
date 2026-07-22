#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::load_date_constraints;
  use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
  use crate::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::partition::timetree::GraphTimetree;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use maplit::btreemap;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::nwk::nwk_read_str;

  const SMALL_TREE_NWK: &str = "((leaf1:1.0,leaf2:1.0)internal1:1.0,leaf3:1.0)root:1.0;";

  fn small_tree_dates() -> DatesMap {
    btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "internal1".to_owned() => Some(DateConstraint::exact(2005.0)),
      "leaf1".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf3".to_owned() => Some(DateConstraint::exact(2012.0)),
    }
  }

  #[test]
  fn test_optimize_skyline_returns_result() -> Result<(), Report> {
    let graph = helpers::create_graph_with_dates(SMALL_TREE_NWK, &small_tree_dates())?;
    let params = SkylineParams {
      n_points: 5,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params)?;

    assert_eq!(result.segment_boundaries.len(), 6);
    assert_eq!(result.tc_values.len(), 5);
    assert!(result.log_likelihood.is_finite());

    Ok(())
  }

  #[test]
  fn test_optimize_skyline_tc_values_positive() -> Result<(), Report> {
    let graph = helpers::create_graph_with_dates(SMALL_TREE_NWK, &small_tree_dates())?;
    let params = SkylineParams {
      n_points: 5,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params)?;

    for &tc in &result.tc_values {
      assert!(
        tc > 0.0 && tc.is_finite(),
        "Tc segment value must be positive and finite, got {tc}"
      );
    }

    Ok(())
  }

  #[test]
  fn test_optimize_skyline_tc_distribution_evaluates() -> Result<(), Report> {
    let graph = helpers::create_graph_with_dates(SMALL_TREE_NWK, &small_tree_dates())?;
    let params = SkylineParams {
      n_points: 5,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params)?;

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

    let graph = helpers::create_graph_with_dates(TREE_NWK, &dates)?;
    let params = SkylineParams {
      n_points: 10,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params)?;

    assert_eq!(result.tc_values.len(), 10);
    assert!(result.log_likelihood.is_finite());

    Ok(())
  }

  #[test]
  fn test_skyline_reported_likelihood_matches_model_evaluation() -> Result<(), Report> {
    // The reported log-likelihood must equal the shared per-edge model cost
    // evaluated on the returned piecewise-constant Tc(t).
    let graph = helpers::create_graph_with_dates(SMALL_TREE_NWK, &small_tree_dates())?;
    let params = SkylineParams {
      n_points: 4,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params)?;
    let expected = compute_coalescent_total_lh(&graph, &result.tc_distribution)?;

    pretty_assert_ulps_eq!(expected, result.log_likelihood, max_ulps = 10);
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
    let graph = helpers::create_graph_with_dates(TREE_NWK, &dates)?;
    let params = SkylineParams {
      n_points: 2,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params)?;
    let expected = compute_coalescent_total_lh(&graph, &result.tc_distribution)?;

    // Oracle: the canonical per-edge Kingman cost assigns m - 1 merger-rate
    // factors to an m-child polytomy (Kingman 1982, doi:10.1016/0304-4149(82)90011-4).
    pretty_assert_ulps_eq!(expected, result.log_likelihood, max_ulps = 10);
    Ok(())
  }

  #[test]
  fn test_skyline_beats_or_matches_constant_tc() -> Result<(), Report> {
    // The regularized skyline optimum should not have lower likelihood than the
    // best constant Tc (a skyline with all segments equal is a feasible point).
    let graph = helpers::create_graph_with_dates(SMALL_TREE_NWK, &small_tree_dates())?;
    let params = SkylineParams {
      n_points: 4,
      stiffness: 1e-6, // near-zero smoothing: skyline free to fit each segment
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params)?;

    let constant_tc = crate::coalescent::optimize_tc::optimize_tc(&graph)?;
    assert!(
      result.log_likelihood >= constant_tc.likelihood - 1e-6,
      "skyline LL {} should be >= constant-Tc LL {}",
      result.log_likelihood,
      constant_tc.likelihood
    );

    Ok(())
  }

  mod helpers {
    use super::*;

    pub fn create_graph_with_dates(tree_nwk: &str, dates: &DatesMap) -> Result<GraphTimetree, Report> {
      let graph = nwk_read_str(tree_nwk)?;
      load_date_constraints(dates, &graph)?;
      Ok(graph)
    }
  }
}
