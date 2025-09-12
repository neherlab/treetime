use crate::commands::clock::clock_graph::ClockGraph;
use crate::commands::clock::clock_regression::ClockOptions;
use crate::commands::clock::find_best_root::find_best_split::{FindRootResult, find_best_split};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::utils::container::get_exactly_one;
use eyre::Report;
use std::sync::Arc;

/// Find the best new root node
///
// Loop over all nodes, pick the one with the lowest chisq, then optimize position along surrounding branches.
pub fn find_best_root(
  graph: &ClockGraph,
  options: &ClockOptions,
  params: &BranchPointOptimizationParams,
) -> Result<FindRootResult, Report> {
  let root = graph.get_exactly_one_root()?;
  let mut best_root_node = Arc::clone(&root);

  // Initialize with the current root
  let root = root.read_arc().payload().read_arc();
  let mut best_chisq = root.total.chisq();
  let mut best_res = FindRootResult {
    edge: None,
    split: 0.0,
    chisq: best_chisq,
    total: root.total.clone(),
  };

  // Find best node
  for n in graph.get_nodes() {
    let tmp_chisq = n.read_arc().payload().read_arc().total.chisq();
    if tmp_chisq < best_chisq {
      best_chisq = tmp_chisq;
      best_root_node = Arc::clone(&n);
    }
  }
  let best_root_node = best_root_node.read_arc();

  // Check if some intermediate place on the parent branch is better
  if !best_root_node.is_root() {
    let inbound = best_root_node.inbound();
    let edge = get_exactly_one(inbound).expect("Not implemented: multiple parent nodes");
    let res = find_best_split(graph, *edge, options, params)?;
    if res.chisq < best_chisq {
      best_chisq = res.chisq;
      best_res = res;
    }
  }

  // Check if some place on a child branch is better
  for e in best_root_node.outbound() {
    let res = find_best_split(graph, *e, options, params)?;
    if res.chisq < best_chisq {
      best_chisq = res.chisq;
      best_res = res;
    }
  }

  Ok(best_res)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::commands::clock::clock_regression::{clock_regression_backward, clock_regression_forward};
  use crate::commands::clock::find_best_root::params::{BrentParams, GoldenSectionParams, GridSearchParams};
  use crate::graph::node::Named;
  use crate::io::nwk::nwk_read_str;
  use crate::o;
  use approx::assert_ulps_eq;
  use maplit::btreemap;

  fn setup_test_graph() -> Result<(ClockGraph, ClockOptions), Report> {
    let dates = btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
    };

    let graph: ClockGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    for n in graph.get_leaves() {
      let name = n.read_arc().payload().read_arc().name().unwrap().as_ref().to_owned();
      n.write_arc().payload().write_arc().date = Some(dates[&name]);
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
        max_iters: 25,
        tolerance: 1e-8,
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
        max_iters: 25,
        tolerance: 1e-8,
      }),
    )?;

    assert_ulps_eq!(best_root.chisq, 0.00025599996471386156, epsilon = 1e-16);

    Ok(())
  }
}
