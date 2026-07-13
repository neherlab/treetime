#[cfg(test)]
mod tests {
  use crate::payload::ancestral::GraphAncestral;
  use crate::reroot::div_stats::DivStats;
  use crate::reroot::div_stats_traversal::compute_div_stats;
  use crate::reroot::params::BranchPointOptimizationParams;
  use crate::reroot::search::find_best_root;
  use crate::reroot::traits::RootStats;
  use crate::reroot::variance::VarianceModel;
  use approx::{assert_abs_diff_eq, assert_ulps_eq};
  use eyre::Report;
  use treetime_io::nwk::nwk_read_str;

  // Root-to-tip distances {0.1, 0.3} => variance (0.1-0.3)^2/4 = 0.01 at the
  // current root (default variance model: leaf var 1, internal var 0).
  #[test]
  fn test_search_root_stats_score_matches_analytical() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.3)root;")?;
    let field = compute_div_stats(&graph, &VarianceModel::default())?;
    assert_ulps_eq!(field.root_stats.count(), 2.0, max_ulps = 4);
    assert_ulps_eq!(field.root_stats.d_sum(), 0.4, max_ulps = 4);
    assert_ulps_eq!(field.root_stats.dsq_sum(), 0.10, max_ulps = 4);
    assert_ulps_eq!(field.root_stats.score(), 0.01, max_ulps = 8);
    Ok(())
  }

  // The unbalanced root is improvable: shifting the root toward the longer branch
  // makes both tips equidistant (variance 0). Brent finds the interior optimum.
  #[test]
  fn test_search_finds_equidistant_root_on_unbalanced_tree() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.3)root;")?;
    let field = compute_div_stats(&graph, &VarianceModel::default())?;
    let best = find_best_root(
      &graph,
      &field.edge_stats,
      &field.root_stats,
      &VarianceModel::default(),
      &BranchPointOptimizationParams::brent(),
    )?;

    assert!(best.edge.is_some(), "expected an improving reroot");
    assert_abs_diff_eq!(best.score, 0.0, epsilon = 1e-10);
    assert!(best.score < field.root_stats.score());
    Ok(())
  }

  // An already-equidistant root cannot be improved: the baseline (edge = None) wins.
  #[test]
  fn test_search_keeps_balanced_root() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.2,B:0.2)root;")?;
    let field = compute_div_stats(&graph, &VarianceModel::default())?;
    let best = find_best_root(
      &graph,
      &field.edge_stats,
      &field.root_stats,
      &VarianceModel::default(),
      &BranchPointOptimizationParams::brent(),
    )?;

    assert!(best.edge.is_none(), "balanced root should not be improved");
    assert_ulps_eq!(best.score, 0.0, max_ulps = 8);
    Ok(())
  }

  // A symmetric star tree is already optimal regardless of arity.
  #[test]
  fn test_search_keeps_balanced_star_root() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.2,B:0.2,C:0.2)root;")?;
    let field = compute_div_stats(&graph, &VarianceModel::default())?;
    let best = find_best_root(
      &graph,
      &field.edge_stats,
      &field.root_stats,
      &VarianceModel::default(),
      &BranchPointOptimizationParams::brent(),
    )?;

    assert!(best.edge.is_none());
    assert_ulps_eq!(field.root_stats.score(), 0.0, max_ulps = 8);
    Ok(())
  }

  // Every edge must receive both directional messages from the traversal.
  #[test]
  fn test_search_traversal_covers_all_edges() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)i:0.3,C:0.4)root;")?;
    let field = compute_div_stats(&graph, &VarianceModel::default())?;
    assert_eq!(field.edge_stats.len(), graph.get_edges().len());
    for (to_parent, to_child) in field.edge_stats.values() {
      assert!(to_parent.count() > 0.0);
      assert!(to_child.count() > 0.0);
    }
    Ok(())
  }

  #[test]
  fn test_search_root_stats_trait_dispatch() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.3)root;")?;
    let field = compute_div_stats(&graph, &VarianceModel::default())?;
    let score = <DivStats as RootStats>::score(&field.root_stats);
    assert_ulps_eq!(score, 0.01, max_ulps = 8);
    Ok(())
  }
}
