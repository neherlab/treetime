#[cfg(test)]
mod tests {
  use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
  use crate::payload::traits::TimetreeNode;
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::inference::forward_pass::propagate_distributions_forward;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::node::TimeConstraint;
  use treetime_io::nwk::nwk_read_str;

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
    graph
      .get_node(leaf_key)
      .expect("leaf A exists")
      .read_arc()
      .payload()
      .write_arc()
      .set_time_distribution(Some(Arc::new(Distribution::point(2013.0, 1.0))));

    propagate_distributions_forward(&graph)?;

    let root_time = graph
      .get_node(root_key)
      .expect("root exists")
      .read_arc()
      .payload()
      .read_arc()
      .time();
    let leaf_time = graph
      .get_node(leaf_key)
      .expect("leaf A exists")
      .read_arc()
      .payload()
      .read_arc()
      .time();

    assert_eq!(None, root_time);
    let leaf_time = leaf_time.expect("leaf A should keep its observed date");
    pretty_assert_ulps_eq!(leaf_time, 2013.0, max_ulps = 4);

    Ok(())
  }
}
