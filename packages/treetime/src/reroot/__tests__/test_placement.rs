#[cfg(test)]
mod tests {
  use crate::reroot::placement::{RootTarget, leaf_keys, require_dated_new_leaves};
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::reroot::StemRemovalInfo;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::assert_error;

  use helpers::{Expected, edge_fixture};

  #[rustfmt::skip]
  #[rstest]
  #[case::zero(                    0.0,                        true,  Expected::Source)]
  #[case::within_epsilon_of_zero(  f64::EPSILON / 2.0,         true,  Expected::Source)]
  #[case::one(                     1.0,                        true,  Expected::Target)]
  #[case::within_epsilon_of_one(   1.0 - f64::EPSILON / 2.0,   true,  Expected::Target)]
  #[case::interior_split(          0.25,                       true,  Expected::Split)]
  #[case::interior_split_upper(    0.75,                       true,  Expected::Split)]
  #[case::below_half_no_split(     0.25,                       false, Expected::Source)]
  #[case::half_no_split(           0.5,                        false, Expected::Target)]
  #[case::above_half_no_split(     0.75,                       false, Expected::Target)]
  #[trace]
  fn test_root_target_on_edge(
    #[case] split: f64,
    #[case] split_edge: bool,
    #[case] expected: Expected,
  ) -> Result<(), Report> {
    let (graph, edge_key, source_key, target_key) = edge_fixture()?;
    let expected = match expected {
      Expected::Source => RootTarget::Node(source_key),
      Expected::Target => RootTarget::Node(target_key),
      Expected::Split => RootTarget::Split { edge_key, split },
    };

    let actual = RootTarget::on_edge(&graph, edge_key, split, split_edge)?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_root_target_after_stem_removal_moves_a_split_of_the_stem_edge_to_the_stem_child() -> Result<(), Report> {
    let (_, edge_key, source_key, target_key) = edge_fixture()?;
    let stem = StemRemovalInfo {
      removed_node_key: source_key,
      removed_edge_key: edge_key,
      new_root_key: target_key,
    };

    let actual = RootTarget::Split { edge_key, split: 0.25 }.after_stem_removal(&stem);

    assert_eq!(RootTarget::Node(target_key), actual);
    Ok(())
  }

  #[test]
  fn test_root_target_after_stem_removal_keeps_targets_off_the_stem_edge() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.3)STEM;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let stem = StemRemovalInfo {
      removed_node_key: find_node_key_by_name(&graph, &names, "STEM").expect("STEM exists"),
      removed_edge_key: find_edge_key(&graph, &names, "STEM", "AB").expect("STEM->AB exists"),
      new_root_key: find_node_key_by_name(&graph, &names, "AB").expect("AB exists"),
    };
    let edge_key = find_edge_key(&graph, &names, "AB", "A").expect("AB->A exists");
    let a_key = find_node_key_by_name(&graph, &names, "A").expect("A exists");

    let split = RootTarget::Split { edge_key, split: 0.25 }.after_stem_removal(&stem);
    let node = RootTarget::Node(a_key).after_stem_removal(&stem);

    assert_eq!(RootTarget::Split { edge_key, split: 0.25 }, split);
    assert_eq!(RootTarget::Node(a_key), node);
    Ok(())
  }

  #[test]
  fn test_require_dated_new_leaves_rejects_a_new_undated_leaf() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.3,C:0.4)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let b_key = find_node_key_by_name(&graph, &names, "B").expect("B exists");
    let mut leaves_before = leaf_keys(&graph);
    leaves_before.remove(&b_key);

    let result = require_dated_new_leaves(&graph, &leaves_before, |_| false, &names);

    assert_error!(
      result,
      "Rerooting turned the internal node 'B' into a leaf, but it has no date and no observed data. Remove the node from the input tree or pass --keep-root."
    );
    Ok(())
  }

  #[test]
  fn test_require_dated_new_leaves_accepts_a_new_dated_leaf() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.3,C:0.4)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let b_key = find_node_key_by_name(&graph, &names, "B").expect("B exists");
    let mut leaves_before = leaf_keys(&graph);
    leaves_before.remove(&b_key);

    require_dated_new_leaves(&graph, &leaves_before, |key| key == b_key, &names)?;
    Ok(())
  }

  #[test]
  fn test_require_dated_new_leaves_ignores_undated_leaves_that_were_leaves_before() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.3,C:0.4)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    require_dated_new_leaves(&graph, &leaf_keys(&graph), |_| false, &names)?;
    Ok(())
  }

  mod helpers {
    use crate::test_utils::{find_edge_key, find_node_key_by_name};
    use eyre::Report;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read_str;

    #[derive(Debug, Clone, Copy)]
    pub(super) enum Expected {
      Source,
      Target,
      Split,
    }

    pub(super) fn edge_fixture() -> Result<(Graph, GraphEdgeKey, GraphNodeKey, GraphNodeKey), Report> {
      let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.3,C:0.4)root;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let edge_key = find_edge_key(&graph, &names, "root", "AB").expect("root->AB exists");
      let source_key = find_node_key_by_name(&graph, &names, "root").expect("root exists");
      let target_key = find_node_key_by_name(&graph, &names, "AB").expect("AB exists");
      Ok((graph, edge_key, source_key, target_key))
    }
  }
}
