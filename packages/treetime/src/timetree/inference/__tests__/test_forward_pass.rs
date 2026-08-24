#[cfg(test)]
mod tests {
  use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
  use crate::payload::traits::TimetreeNode;
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::inference::forward_pass::{propagate_distributions_forward, set_likely_time};
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::edge::BranchDistribution;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::{GraphNodeKey, TimeConstraint};
  use treetime_io::nwk::nwk_read_str;

  type TestGraph = Graph<NodeTimetree, EdgeTimetree, ()>;

  #[test]
  fn test_forward_pass_set_likely_time_empty_distribution_returns_none() {
    let mut node = node_with_distribution(Some(Distribution::empty()));
    assert_eq!(None, set_likely_time(&mut node, None));
    assert_eq!(None, node.time());
  }

  #[test]
  fn test_forward_pass_set_likely_time_missing_distribution_returns_none() {
    let mut node = node_with_distribution(None);
    assert_eq!(None, set_likely_time(&mut node, None));
    assert_eq!(None, node.time());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_parent(      None,      5.0)]
  #[case::parent_earlier( Some(3.0), 5.0)]
  #[case::parent_clamps(  Some(8.0), 8.0)]
  #[trace]
  fn test_forward_pass_set_likely_time_uses_distribution_peak(
    #[case] parent_time: Option<f64>,
    #[case] expected: f64,
  ) {
    let mut node = node_with_distribution(Some(Distribution::point(5.0, 1.0)));
    let assigned = set_likely_time(&mut node, parent_time).expect("a time should be assigned");
    pretty_assert_ulps_eq!(assigned, expected, max_ulps = 4);
    let committed = node.time().expect("node time should be committed");
    pretty_assert_ulps_eq!(committed, expected, max_ulps = 4);
  }

  /// The forward pass must not assign a date to an internal node whose time distribution is
  /// empty (irreconcilable constraints). It leaves the node undated and warns, rather than
  /// silently emitting a structurally valid tree with missing internal dates. Leaves keep their
  /// observed dates.
  #[test]
  fn test_forward_pass_leaves_internal_node_with_empty_distribution_undated() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("(A:2.5)root;")?;

    let root_key = find_node_key_by_name(&graph, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");

    // Internal root carries an empty (irreconcilable) time distribution; the leaf carries a date.
    graph
      .get_node(root_key)
      .expect("root exists")
      .read_arc()
      .payload()
      .write_arc()
      .set_time_distribution(Some(Arc::new(Distribution::empty())));
    set_date(&graph, leaf_key, Distribution::point(2013.0, 0.0));

    propagate_distributions_forward(&graph)?;

    let root_time = graph
      .get_node(root_key)
      .expect("root exists")
      .read_arc()
      .payload()
      .read_arc()
      .time();

    assert_eq!(None, root_time);
    let leaf_time = node_time(&graph, leaf_key).expect("leaf A should keep its observed date");
    pretty_assert_ulps_eq!(leaf_time, 2013.0, max_ulps = 4);

    Ok(())
  }

  /// A leaf whose date is uncertain is refined from its parent exactly as an internal node is: the
  /// date range says where the leaf can be, and the message coming down the tree says where within
  /// that range it most likely is. Here the root sits at 2009 and the branch to the leaf takes one
  /// year, so the leaf lands at 2010 rather than at 2011.5, the midpoint of the range it was given.
  #[test]
  fn test_forward_pass_refines_uncertain_leaf_date_from_parent() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("(A:1.0)root;")?;
    let root_key = find_node_key_by_name(&graph, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");

    set_time_distribution(&graph, root_key, Distribution::point(2009.0, 0.0));
    set_date(&graph, leaf_key, Distribution::range((2009.5, 2013.5), 0.0));
    set_branch_length_distribution(&graph, leaf_key, 1.0);

    propagate_distributions_forward(&graph)?;

    let leaf_time = node_time(&graph, leaf_key).expect("leaf A should be dated");
    pretty_assert_ulps_eq!(leaf_time, 2010.0, max_ulps = 4);

    // The refinement is on the distribution, not just on the point estimate: what the leaf carries
    // afterwards is the posterior the confidence interval is read from.
    let leaf_dist = leaf_time_distribution(&graph, leaf_key).expect("leaf A should have a distribution");
    assert_eq!(Distribution::point(2010.0, 0.0), leaf_dist);

    Ok(())
  }

  /// A leaf given an exact date is not refined and not clamped: the date is an observation, and a
  /// leaf dated before its parent is a conflict between that observation and the fitted clock that
  /// `commit_clock_branch_lengths` reports. Clamping it to the parent would hide the conflict.
  #[test]
  fn test_forward_pass_keeps_exact_leaf_date_earlier_than_its_parent() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("(A:1.0)root;")?;
    let root_key = find_node_key_by_name(&graph, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");

    set_time_distribution(&graph, root_key, Distribution::point(2009.0, 0.0));
    set_date(&graph, leaf_key, Distribution::point(2008.0, 0.0));
    set_branch_length_distribution(&graph, leaf_key, 1.0);

    propagate_distributions_forward(&graph)?;

    let leaf_time = node_time(&graph, leaf_key).expect("leaf A should keep its observed date");
    pretty_assert_ulps_eq!(leaf_time, 2008.0, max_ulps = 4);

    let leaf_dist = leaf_time_distribution(&graph, leaf_key).expect("leaf A should have a distribution");
    assert_eq!(Distribution::point(2008.0, 0.0), leaf_dist);

    Ok(())
  }

  /// A leaf whose date is uncertain respects the parent constraint like any inferred node. Here the
  /// branch carries no length distribution, so there is nothing to refine the range with, and the
  /// peak of the range alone would put the leaf 3 years before its parent.
  #[test]
  fn test_forward_pass_clamps_uncertain_leaf_date_to_parent_time() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("(A:1.0)root;")?;
    let root_key = find_node_key_by_name(&graph, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");

    set_time_distribution(&graph, root_key, Distribution::point(2009.0, 0.0));
    set_date(&graph, leaf_key, Distribution::range((2005.0, 2007.0), 0.0));

    propagate_distributions_forward(&graph)?;

    let leaf_time = node_time(&graph, leaf_key).expect("leaf A should be dated");
    pretty_assert_ulps_eq!(leaf_time, 2009.0, max_ulps = 4);

    Ok(())
  }

  /// A date the tree contradicts outright is kept, not refined away. The parent sits at 2009 and
  /// the branch takes a year, so the message from the parent is a point at 2010 that the leaf's
  /// range [2005, 2007] rules out entirely and their product is empty. Refining onto that would
  /// leave the leaf undated, which is worse than the date the input gave it.
  #[test]
  fn test_forward_pass_keeps_uncertain_leaf_date_the_tree_contradicts() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("(A:1.0)root;")?;
    let root_key = find_node_key_by_name(&graph, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");

    set_time_distribution(&graph, root_key, Distribution::point(2009.0, 0.0));
    let given = Distribution::range((2005.0, 2007.0), 0.0);
    set_date(&graph, leaf_key, given.clone());
    set_branch_length_distribution(&graph, leaf_key, 1.0);

    propagate_distributions_forward(&graph)?;

    let leaf_dist = leaf_time_distribution(&graph, leaf_key).expect("leaf A should keep its given date");
    assert_eq!(given, leaf_dist);

    // Still projected onto the parent, as any node whose time the tree infers is.
    let leaf_time = node_time(&graph, leaf_key).expect("leaf A should be dated");
    pretty_assert_ulps_eq!(leaf_time, 2009.0, max_ulps = 4);

    Ok(())
  }

  /// Cutting the parent's grid down to the window the branch can carry into the leaf's date range
  /// is an optimization, so it must not move the answer. The parent's posterior spans a century at
  /// a resolution that makes it far wider than the leaf's month-wide range.
  #[test]
  fn test_forward_pass_refinement_is_unchanged_by_the_reachable_window() -> Result<(), Report> {
    let wide_parent = {
      let t = Array1::linspace(1950.0, 2050.0, 4001);
      // A Gaussian peaked at 2009 stored on the neg-log axis: the ordinate is `-ln p = 0.5 z^2`, a
      // parabola whose minimum (the peak) sits at 2009, so the leaf's own range decides where in it
      // the leaf lands.
      let y = t.mapv(|t: f64| 0.5 * ((t - 2009.0) / 3.0).powi(2));
      Distribution::function(t, y)?
    };

    let refine = |parent: Distribution<NegLog>| -> Result<f64, Report> {
      let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("(A:1.0)root;")?;
      let root_key = find_node_key_by_name(&graph, "root").expect("root not found");
      let leaf_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");
      set_time_distribution(&graph, root_key, parent);
      set_date(&graph, leaf_key, Distribution::range((2010.5, 2010.6), 0.0));
      set_branch_length_distribution(&graph, leaf_key, 1.0);
      propagate_distributions_forward(&graph)?;
      node_time(&graph, leaf_key).ok_or_else(|| eyre::eyre!("leaf A should be dated"))
    };

    // The same parent, once as given and once already narrowed to a window the restriction cannot
    // cut further: both must date the leaf identically.
    let windowed = match &wide_parent {
      Distribution::Function(f) => Distribution::Function(f.resample_range_dx((2009.4, 2009.7), f.dx())?),
      other => other.clone(),
    };

    pretty_assert_ulps_eq!(refine(wide_parent)?, refine(windowed)?, max_ulps = 4);

    Ok(())
  }

  /// A date range narrower than the parent's own grid resolution still has to be refinable: the
  /// window cut around it cannot be gridded at that resolution, so it is widened rather than used
  /// as is. The parent here resolves a year, the date range a single day.
  #[test]
  fn test_forward_pass_refines_a_date_range_narrower_than_the_parent_grid() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("(A:1.0)root;")?;
    let root_key = find_node_key_by_name(&graph, "root").expect("root not found");
    let leaf_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");

    let coarse_parent = {
      let t = Array1::linspace(1900.0, 2020.0, 121); // one grid point per year
      // Gaussian peaked at 2009 on the neg-log axis (`-ln p = 0.5 z^2`); see the wide-parent case.
      let y = t.mapv(|t: f64| 0.5 * ((t - 2009.0) / 5.0).powi(2));
      Distribution::function(t, y)?
    };
    set_time_distribution(&graph, root_key, coarse_parent);
    set_date(&graph, leaf_key, Distribution::range((2010.500, 2010.503), 0.0));
    set_branch_length_distribution(&graph, leaf_key, 1.0);

    propagate_distributions_forward(&graph)?;

    let leaf_time = node_time(&graph, leaf_key).expect("leaf A should be dated");
    assert!(
      (2010.500..=2010.503).contains(&leaf_time),
      "leaf A should be dated within the day it was given, got {leaf_time}"
    );

    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn node_with_distribution(distribution: Option<Distribution<NegLog>>) -> NodeTimetree {
      NodeTimetree {
        time_distribution: distribution.map(Arc::new),
        ..NodeTimetree::default()
      }
    }

    /// Give a node a date, as loading date constraints from the input does: the fixed constraint
    /// and the time distribution it seeds.
    pub(super) fn set_date(graph: &TestGraph, key: GraphNodeKey, dist: Distribution<NegLog>) {
      let dist = Arc::new(dist);
      let node = graph.get_node(key).expect("node exists");
      let mut payload = node.read_arc().payload().write_arc();
      payload.set_date_constraint(Some(Arc::clone(&dist)));
      payload.set_time_distribution(Some(dist));
    }

    pub(super) fn set_time_distribution(graph: &TestGraph, key: GraphNodeKey, dist: Distribution<NegLog>) {
      graph
        .get_node(key)
        .expect("node exists")
        .read_arc()
        .payload()
        .write_arc()
        .set_time_distribution(Some(Arc::new(dist)));
    }

    pub(super) fn set_branch_length_distribution(graph: &TestGraph, target_key: GraphNodeKey, branch_length: f64) {
      for edge in graph.get_edges() {
        let edge = edge.read_arc();
        if edge.target() == target_key {
          edge
            .payload()
            .write_arc()
            .set_branch_length_distribution(Some(Arc::new(Distribution::point(branch_length, 0.0))));
        }
      }
    }

    pub(super) fn node_time(graph: &TestGraph, key: GraphNodeKey) -> Option<f64> {
      graph
        .get_node(key)
        .expect("node exists")
        .read_arc()
        .payload()
        .read_arc()
        .time()
    }

    pub(super) fn leaf_time_distribution(graph: &TestGraph, key: GraphNodeKey) -> Option<Distribution<NegLog>> {
      graph
        .get_node(key)
        .expect("node exists")
        .read_arc()
        .payload()
        .read_arc()
        .time_distribution()
        .as_ref()
        .map(|dist| dist.as_ref().clone())
    }
  }

  use helpers::*;
}
