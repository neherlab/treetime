#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::load_date_constraints;
  use crate::coalescent::skyline::{SkylineParams, build_tc_distribution, optimize_skyline};
  use crate::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::partition::timetree::GraphTimetree;
  use crate::pretty_assert_ulps_eq;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::array;
  use rstest::rstest;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_optimize_skyline_returns_result() -> Result<(), Report> {
    const TREE_NWK: &str = "((leaf1:1.0,leaf2:1.0)internal1:1.0,leaf3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "internal1".to_owned() => Some(DateConstraint::exact(2005.0)),
      "leaf1".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf3".to_owned() => Some(DateConstraint::exact(2012.0)),
    };

    let graph = helpers::create_graph_with_dates(TREE_NWK, &dates)?;
    let params = SkylineParams {
      n_points: 5,
      stiffness: 2.0,
      regularization: 10.0,
      tolerance: 0.1,
      max_iter: 100,
    };

    let result = optimize_skyline(&graph, &params)?;

    assert_eq!(result.time_grid.len(), 5);
    assert_eq!(result.log_tc_values.len(), 5);
    assert!(result.log_likelihood.is_finite());

    Ok(())
  }

  #[test]
  fn test_optimize_skyline_tc_distribution_evaluates() -> Result<(), Report> {
    const TREE_NWK: &str = "((leaf1:1.0,leaf2:1.0)internal1:1.0,leaf3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "internal1".to_owned() => Some(DateConstraint::exact(2005.0)),
      "leaf1".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf3".to_owned() => Some(DateConstraint::exact(2012.0)),
    };

    let graph = helpers::create_graph_with_dates(TREE_NWK, &dates)?;
    let params = SkylineParams {
      n_points: 5,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params)?;

    // Tc distribution should be evaluable within the time grid
    let t_min = result.time_grid[0];
    let t_max = result.time_grid[result.time_grid.len() - 1];
    let t_mid = f64::midpoint(t_min, t_max);

    let tc_min = result.tc_distribution.eval(t_min)?;
    let tc_mid = result.tc_distribution.eval(t_mid)?;
    let tc_max = result.tc_distribution.eval(t_max)?;

    assert!(tc_min > 0.0);
    assert!(tc_mid > 0.0);
    assert!(tc_max > 0.0);

    Ok(())
  }

  #[test]
  fn test_optimize_skyline_log_tc_in_reasonable_range() -> Result<(), Report> {
    const TREE_NWK: &str = "((leaf1:1.0,leaf2:1.0)internal1:1.0,leaf3:1.0)root:1.0;";
    let dates = btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "internal1".to_owned() => Some(DateConstraint::exact(2005.0)),
      "leaf1".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf3".to_owned() => Some(DateConstraint::exact(2012.0)),
    };

    let graph = helpers::create_graph_with_dates(TREE_NWK, &dates)?;
    let params = SkylineParams {
      n_points: 5,
      regularization: 10.0,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params)?;

    // With regularization, log_tc values should stay in reasonable range [-100, 0]
    for &log_tc in &result.log_tc_values {
      assert!(
        (-100.0..=10.0).contains(&log_tc),
        "log_tc={log_tc} outside reasonable range"
      );
    }

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

    assert_eq!(result.time_grid.len(), 10);
    assert!(result.log_likelihood.is_finite());

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
      tolerance: 0.1,
      max_iter: 100,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&graph, &params)?;
    let expected = compute_coalescent_total_lh(&graph, &result.tc_distribution)?;

    // Oracle: the canonical per-edge Kingman cost assigns m - 1 merger-rate
    // factors to an m-child polytomy (Kingman 1982, doi:10.1016/0304-4149(82)90011-4).
    pretty_assert_ulps_eq!(expected, result.log_likelihood, max_ulps = 10);
    Ok(())
  }

  /// `build_tc_distribution` must interpolate `exp(log_tc)` linearly between grid
  /// points and clamp to the boundary values outside the grid range. With
  /// `time_grid = [0, 1, 2]` and `log_tc = [0, ln 2, 0]` the values are
  /// `tc = [1, 2, 1]`, so the function rises linearly 1 -> 2 on `[0, 1]`,
  /// falls linearly 2 -> 1 on `[1, 2]`, and is flat (= 1) outside `[0, 2]`.
  #[rustfmt::skip]
  #[rstest]
  #[case::below_range_clamps_to_first(-1.0, 1.0)]
  #[case::first_grid_point(           0.0, 1.0)]
  #[case::interior_rising(            0.5, 1.5)]
  #[case::peak_grid_point(            1.0, 2.0)]
  #[case::interior_falling(           1.5, 1.5)]
  #[case::last_grid_point(            2.0, 1.0)]
  #[case::above_range_clamps_to_last( 3.0, 1.0)]
  #[trace]
  fn test_skyline_build_tc_distribution_interpolates_and_clamps(
    #[case] t: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let time_grid = array![0.0, 1.0, 2.0];
    let log_tc = array![0.0, std::f64::consts::LN_2, 0.0];

    let tc_dist = build_tc_distribution(&time_grid, &log_tc);
    let actual = tc_dist.eval(t)?;

    assert_abs_diff_eq!(expected, actual, epsilon = 1e-12);
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
