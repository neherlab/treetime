#[cfg(test)]
mod tests {
  use crate::clock::clock_graph::GraphClock;
  use crate::clock::clock_regression::{ClockParams, clock_regression_backward, clock_regression_forward};
  use crate::clock::find_best_root::find_best_root::find_best_root;
  use crate::clock::find_best_root::find_best_split::FindRootResult;
  use crate::clock::find_best_root::params::{
    BranchPointOptimizationParams, BrentParams, GoldenSectionParams, GridSearchParams, RootObjective,
  };
  use crate::o;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use maplit::btreemap;
  use std::collections::BTreeMap;
  use treetime_graph::node::Named;
  use treetime_io::nwk::nwk_read_str;

  fn setup_graph_with_dates(dates: &BTreeMap<String, f64>) -> Result<(GraphClock, ClockParams), Report> {
    let graph: GraphClock = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    for n in graph.get_leaves() {
      let name = n.read_arc().payload().read_arc().name().unwrap().as_ref().to_owned();
      n.write_arc().payload().write_arc().time = Some(dates[&name]);
    }

    let options = ClockParams::default();
    clock_regression_backward(&graph, &options, None)?;
    clock_regression_forward(&graph, &options, None)?;

    Ok((graph, options))
  }

  fn setup_test_graph() -> Result<(GraphClock, ClockParams), Report> {
    let dates = btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
    };
    setup_graph_with_dates(&dates)
  }

  fn get_edge_node_names(graph: &GraphClock, result: &FindRootResult) -> (String, String) {
    let edge_key = result.edge.expect("result should have an edge");
    let edge = graph.get_edge(edge_key).expect("edge should exist");
    let edge = edge.read_arc();
    let source = graph.get_node(edge.source()).expect("source should exist");
    let target = graph.get_node(edge.target()).expect("target should exist");
    let source_name = source
      .read_arc()
      .payload()
      .read_arc()
      .name()
      .map_or_else(|| "unnamed".to_owned(), |n| n.as_ref().to_owned());
    let target_name = target
      .read_arc()
      .payload()
      .read_arc()
      .name()
      .map_or_else(|| "unnamed".to_owned(), |n| n.as_ref().to_owned());
    (source_name, target_name)
  }

  #[test]
  fn test_find_best_root_grid() -> Result<(), Report> {
    let (graph, options) = setup_test_graph()?;

    let best_root = find_best_root(
      &graph,
      &options,
      &BranchPointOptimizationParams::grid(),
      true,
      RootObjective::EstimatedRate,
    )?;

    // Verify chisq value
    pretty_assert_ulps_eq!(best_root.chisq, 0.0002610661988682317, max_ulps = 4);

    // Verify split is in valid range
    assert!(
      best_root.split >= 0.0 && best_root.split <= 1.0,
      "split should be in [0, 1]"
    );

    // Verify edge position - best root is on edge from root to CD
    let (source, target) = get_edge_node_names(&graph, &best_root);
    assert_eq!(source, "root");
    assert_eq!(target, "CD");

    Ok(())
  }

  #[test]
  fn test_find_best_root_grid_with_params() -> Result<(), Report> {
    let (graph, options) = setup_test_graph()?;

    let best_root = find_best_root(
      &graph,
      &options,
      &BranchPointOptimizationParams::grid_with(GridSearchParams { n_points: 51 }),
      true,
      RootObjective::EstimatedRate,
    )?;

    // Verify chisq value
    pretty_assert_ulps_eq!(best_root.chisq, 0.0002560258129903322, max_ulps = 4);

    // Verify split is in valid range
    assert!(
      best_root.split >= 0.0 && best_root.split <= 1.0,
      "split should be in [0, 1]"
    );

    // Verify edge position - best root is on edge from root to CD
    let (source, target) = get_edge_node_names(&graph, &best_root);
    assert_eq!(source, "root");
    assert_eq!(target, "CD");

    Ok(())
  }

  #[test]
  fn test_find_best_root_grid_scores_fixed_rate_objective() -> Result<(), Report> {
    let (graph, options) = setup_test_graph()?;

    let best_root = find_best_root(
      &graph,
      &options,
      &BranchPointOptimizationParams::grid(),
      true,
      RootObjective::FixedRate(0.0),
    )?;

    let expected_chisq = best_root.clock_set.chisq_fixed_rate(0.0);
    pretty_assert_ulps_eq!(best_root.chisq, expected_chisq, max_ulps = 4);
    assert!(
      (best_root.chisq - best_root.clock_set.chisq()).abs() > 1e-10,
      "test must distinguish fixed-rate and estimated-rate objectives"
    );

    Ok(())
  }

  #[test]
  fn test_find_best_root_brent() -> Result<(), Report> {
    let (graph, options) = setup_test_graph()?;

    let best_root = find_best_root(
      &graph,
      &options,
      &BranchPointOptimizationParams::brent(),
      true,
      RootObjective::EstimatedRate,
    )?;

    // Verify chisq value
    pretty_assert_ulps_eq!(best_root.chisq, 0.00025599996471448085, max_ulps = 4);

    // Verify split is in valid range
    assert!(
      best_root.split >= 0.0 && best_root.split <= 1.0,
      "split should be in [0, 1]"
    );

    // Verify edge position - best root is on edge from root to CD
    let (source, target) = get_edge_node_names(&graph, &best_root);
    assert_eq!(source, "root");
    assert_eq!(target, "CD");

    Ok(())
  }

  #[test]
  fn test_find_best_root_brent_with_params() -> Result<(), Report> {
    let (graph, options) = setup_test_graph()?;

    let best_root = find_best_root(
      &graph,
      &options,
      &BranchPointOptimizationParams::brent_with(BrentParams {
        brent_max_iters: 25,
        brent_tolerance: 1e-8,
      }),
      true,
      RootObjective::EstimatedRate,
    )?;

    // Verify chisq value
    pretty_assert_ulps_eq!(best_root.chisq, 0.00025599996471448085, max_ulps = 4);

    // Verify split is in valid range
    assert!(
      best_root.split >= 0.0 && best_root.split <= 1.0,
      "split should be in [0, 1]"
    );

    // Verify edge position - best root is on edge from root to CD
    let (source, target) = get_edge_node_names(&graph, &best_root);
    assert_eq!(source, "root");
    assert_eq!(target, "CD");

    Ok(())
  }

  #[test]
  fn test_find_best_root_golden_section() -> Result<(), Report> {
    let (graph, options) = setup_test_graph()?;

    let best_root = find_best_root(
      &graph,
      &options,
      &BranchPointOptimizationParams::golden_section(),
      true,
      RootObjective::EstimatedRate,
    )?;

    // Verify chisq value
    pretty_assert_ulps_eq!(best_root.chisq, 0.00025599996471244515, max_ulps = 4);

    // Verify split is in valid range
    assert!(
      best_root.split >= 0.0 && best_root.split <= 1.0,
      "split should be in [0, 1]"
    );

    // Verify edge position - best root is on edge from root to CD
    let (source, target) = get_edge_node_names(&graph, &best_root);
    assert_eq!(source, "root");
    assert_eq!(target, "CD");

    Ok(())
  }

  #[test]
  fn test_find_best_root_golden_section_with_params() -> Result<(), Report> {
    let (graph, options) = setup_test_graph()?;

    let best_root = find_best_root(
      &graph,
      &options,
      &BranchPointOptimizationParams::golden_section_with(GoldenSectionParams {
        golden_max_iters: 25,
        golden_tolerance: 1e-8,
      }),
      true,
      RootObjective::EstimatedRate,
    )?;

    // Verify chisq value
    pretty_assert_ulps_eq!(best_root.chisq, 0.00025599996471386156, max_ulps = 4);

    // Verify split is in valid range
    assert!(
      best_root.split >= 0.0 && best_root.split <= 1.0,
      "split should be in [0, 1]"
    );

    // Verify edge position - best root is on edge from root to CD
    let (source, target) = get_edge_node_names(&graph, &best_root);
    assert_eq!(source, "root");
    assert_eq!(target, "CD");

    Ok(())
  }

  #[test]
  fn test_optimization_methods_improve_on_grid() -> Result<(), Report> {
    let (graph, options) = setup_test_graph()?;

    // Run all three methods
    let grid_result = find_best_root(
      &graph,
      &options,
      &BranchPointOptimizationParams::grid(),
      true,
      RootObjective::EstimatedRate,
    )?;
    let brent_result = find_best_root(
      &graph,
      &options,
      &BranchPointOptimizationParams::brent(),
      true,
      RootObjective::EstimatedRate,
    )?;
    let golden_result = find_best_root(
      &graph,
      &options,
      &BranchPointOptimizationParams::golden_section(),
      true,
      RootObjective::EstimatedRate,
    )?;

    // Brent and golden section should find better or equal chisq than default grid
    // (lower chisq is better)
    assert!(
      brent_result.chisq <= grid_result.chisq,
      "Brent ({:.6e}) should achieve <= chisq than grid ({:.6e})",
      brent_result.chisq,
      grid_result.chisq
    );
    assert!(
      golden_result.chisq <= grid_result.chisq,
      "Golden section ({:.6e}) should achieve <= chisq than grid ({:.6e})",
      golden_result.chisq,
      grid_result.chisq
    );

    // All methods should find the same edge
    assert_eq!(
      grid_result.edge, brent_result.edge,
      "Brent should find same edge as grid"
    );
    assert_eq!(
      grid_result.edge, golden_result.edge,
      "Golden section should find same edge as grid"
    );

    Ok(())
  }

  /// Dates inversely correlated with divergence: negative clock rate at all root positions.
  /// Root-to-tip: A=0.2, B=0.3, C=0.25, D=0.17
  fn setup_negative_rate_graph() -> Result<(GraphClock, ClockParams), Report> {
    let dates = btreemap! {
      o!("A") => 2017.0,
      o!("B") => 2005.0,
      o!("C") => 2010.0,
      o!("D") => 2022.0,
    };
    setup_graph_with_dates(&dates)
  }

  #[test]
  fn test_find_best_root_force_positive_true_rejects_negative_rate() -> Result<(), Report> {
    let (graph, options) = setup_negative_rate_graph()?;

    let result = find_best_root(
      &graph,
      &options,
      &BranchPointOptimizationParams::grid(),
      true,
      RootObjective::EstimatedRate,
    );

    assert!(
      result.is_err(),
      "force_positive=true should reject all-negative-rate graph"
    );
    let err_msg = result.unwrap_err().to_string();
    assert!(
      err_msg.contains("Clock rate is negative"),
      "Error message should mention negative rate, got: {err_msg}"
    );

    Ok(())
  }

  #[test]
  fn test_find_best_root_force_positive_false_accepts_negative_rate() -> Result<(), Report> {
    let (graph, options) = setup_negative_rate_graph()?;

    let best_root = find_best_root(
      &graph,
      &options,
      &BranchPointOptimizationParams::grid(),
      false,
      RootObjective::EstimatedRate,
    )?;

    // Confirm the fixture produces a negative-rate scenario
    let det = best_root.clock_set.determinant();
    assert!(det > 0.0, "determinant should be positive");
    let rate = best_root.clock_set.clock_rate(det);
    assert!(
      rate < 0.0,
      "rate should be negative for this test graph, got {rate:.6e}"
    );

    // Non-negative chi-squared from valid weighted least squares
    assert!(best_root.chisq >= 0.0, "chisq should be non-negative");
    assert!(best_root.chisq.is_finite(), "chisq should be finite");

    assert!(
      best_root.split >= 0.0 && best_root.split <= 1.0,
      "split should be in [0, 1]"
    );

    Ok(())
  }
}
