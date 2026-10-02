#[cfg(test)]
mod tests {
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::divergence::final_divergences;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  const TREE_WITH_UNNAMED_INTERNAL: &str = "((A:0.1,B:0.2)AB:0.3,C:0.4)root;";

  #[test]
  fn test_final_divergences_without_filter_unnamed_node_gets_zero_known_issue_m_timetree_unnamed_internal_node_divergence_stale()
  -> Result<(), Report> {
    let parsed = nwk_read_str(TREE_WITH_UNNAMED_INTERNAL)?;
    let names = parsed.names();
    let graph = parsed.graph;
    let (names, unnamed) = helpers::names_with_unnamed_internal(&graph, names);
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).unwrap();

    let actual = final_divergences(&graph, &parsed.branch_lengths, &names, None)?;

    let expected = btreemap! {
      key("root") => 0.0,
      unnamed     => 0.0,
      key("A")    => 0.4,
      key("B")    => 0.5,
      key("C")    => 0.4,
    };
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_final_divergences_unnamed_node_keeps_stale_filter_value_known_issue_m_timetree_unnamed_internal_node_divergence_stale()
  -> Result<(), Report> {
    let parsed = nwk_read_str(TREE_WITH_UNNAMED_INTERNAL)?;
    let names = parsed.names();
    let graph = parsed.graph;
    let (names, unnamed) = helpers::names_with_unnamed_internal(&graph, names);
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).unwrap();
    let absent = GraphNodeKey(graph.get_nodes().count() + 10);
    let filter_divergences = btreemap! {
      unnamed  => 0.7,
      key("A") => 9.0,
      absent   => 1.0,
    };

    let actual = final_divergences(&graph, &parsed.branch_lengths, &names, Some(&filter_divergences))?;

    let expected = btreemap! {
      key("root") => 0.0,
      unnamed     => 0.7,
      key("A")    => 0.4,
      key("B")    => 0.5,
      key("C")    => 0.4,
    };
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_final_divergences_keys_equal_the_graph_nodes() -> Result<(), Report> {
    let parsed = nwk_read_str(TREE_WITH_UNNAMED_INTERNAL)?;
    let names = parsed.names();
    let graph = parsed.graph;
    let (names, _) = helpers::names_with_unnamed_internal(&graph, names);
    let filter_divergences = btreemap! { GraphNodeKey(100) => 1.0 };

    let actual = final_divergences(&graph, &parsed.branch_lengths, &names, Some(&filter_divergences))?;

    let expected = graph.get_nodes().map(|node| node.key()).collect::<Vec<_>>();
    assert_eq!(expected, actual.into_keys().collect::<Vec<_>>());
    Ok(())
  }

  #[test]
  fn test_final_divergences_nodes_sharing_a_name_get_their_own_divergence() -> Result<(), Report> {
    let parsed = nwk_read_str("((A:0.1,A:0.2)X:0.1,B:0.3)root;")?;
    let names = parsed.names();
    let graph = parsed.graph;
    let duplicates = graph
      .get_leaves()
      .filter(|leaf| names[&leaf.key()].as_deref() == Some("A"))
      .map(|leaf| leaf.key())
      .collect::<Vec<_>>();
    assert_eq!(2, duplicates.len());

    let actual = final_divergences(&graph, &parsed.branch_lengths, &names, None)?;

    let mut duplicate_divergences = duplicates.iter().map(|key| actual[key]).collect::<Vec<_>>();
    duplicate_divergences.sort_by(f64::total_cmp);
    assert_abs_diff_eq!(0.2, duplicate_divergences[0], epsilon = 1e-12);
    assert_abs_diff_eq!(0.3, duplicate_divergences[1], epsilon = 1e-12);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn names_with_unnamed_internal(
      graph: &Graph,
      mut names: BTreeMap<GraphNodeKey, Option<String>>,
    ) -> (BTreeMap<GraphNodeKey, Option<String>>, GraphNodeKey) {
      let internal = graph
        .get_nodes()
        .find(|node| !node.is_leaf() && !node.is_root())
        .map(|node| node.key())
        .unwrap();
      names.insert(internal, None);
      (names, internal)
    }
  }
}
