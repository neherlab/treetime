use super::*;
  use crate::commands::clock::clock_graph::GraphClock;
  use crate::commands::clock::clock_regression::{clock_regression_backward, clock_regression_forward};
  use crate::commands::clock::find_best_root::params::{BrentParams, GoldenSectionParams, GridSearchParams};
  use crate::graph::node::Named;
  use crate::io::nwk::nwk_read_str;
  use crate::o;
  use approx::assert_ulps_eq;
  use maplit::btreemap;

  fn setup_test_graph() -> Result<(GraphClock, ClockOptions), Report> {
    let dates = btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
    };

    let graph: GraphClock = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    for n in graph.get_leaves() {
      let name = n.read_arc().payload().read_arc().name().unwrap().as_ref().to_owned();
      n.write_arc().payload().write_arc().time = Some(dates[&name]);
    }

    let options = ClockOptions::default();

    clock_regression_backward(&graph, &options);
    clock_regression_forward(&graph, &options);

    Ok((graph, options))
  }

  #[test]
  fn test_find_best_root_grid() -> Result<(), Report> {
    let (graph, options) = setup_test_graph()?;

    let best_root = find_best_root(&graph, &options, &BranchPointOptimizationParams::grid())?;
    assert_ulps_eq!(best_root.chisq, 0.0002610661988682317, epsilon = 1e-9);

    Ok(())
  }

  #[test]
  fn test_find_best_root_grid_with_params() -> Result<(), Report> {
    let (graph, options) = setup_test_graph()?;

    let best_root = find_best_root(
      &graph,
      &options,
      &BranchPointOptimizationParams::grid_with(GridSearchParams { n_points: 51 }),
    )?;

    assert_ulps_eq!(best_root.chisq, 0.0002560258129903322, epsilon = 1e-12);

    Ok(())
  }

  #[test]
  fn test_find_best_root_brent() -> Result<(), Report> {
    let (graph, options) = setup_test_graph()?;

    let best_root = find_best_root(&graph, &options, &BranchPointOptimizationParams::brent())?;
    assert_ulps_eq!(best_root.chisq, 0.00025599996471448085, epsilon = 1e-12);

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
    )?;

    assert_ulps_eq!(best_root.chisq, 0.00025599996471448085, epsilon = 1e-16);

    Ok(())
  }

  #[test]
  fn test_find_best_root_golden_section() -> Result<(), Report> {
    let (graph, options) = setup_test_graph()?;

    let best_root = find_best_root(&graph, &options, &BranchPointOptimizationParams::golden_section())?;
    assert_ulps_eq!(best_root.chisq, 0.00025599996471244515, epsilon = 1e-12);

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
    )?;

    assert_ulps_eq!(best_root.chisq, 0.00025599996471386156, epsilon = 1e-16);

    Ok(())
  }
