#[cfg(test)]
mod tests {
  use super::super::test_graph_support::tests::NamedGraph;
  use crate::edge::GraphEdgeKey;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  use self::helpers::{EdgeNames, edge_keys};

  #[rustfmt::skip]
  #[rstest]
  #[case::chain(
    (vec!["root", "internal", "A"],                     vec![("root", "internal"), ("internal", "A")]),
    ("root", "internal"),
    (vec![("internal", "A")],                           vec![("internal", "A")],                 vec!["A"]),
  )]
  #[case::sibling_leaf(
    (vec!["root", "internal", "B", "A"],                vec![("root", "internal"), ("root", "B"), ("internal", "A")]),
    ("root", "internal"),
    (vec![("root", "B"), ("internal", "A")],            vec![("internal", "A")],                 vec!["B", "A"]),
  )]
  #[case::sibling_subtree(
    (vec!["root", "left", "right", "A", "B", "C", "D"], vec![("root", "left"), ("root", "right"), ("left", "A"), ("left", "B"), ("right", "C"), ("right", "D")]),
    ("root", "left"),
    (vec![("root", "right"), ("left", "A"), ("left", "B")], vec![("left", "A"), ("left", "B")],    vec!["right", "A", "B"]),
  )]
  #[case::leaf_edge(
    (vec!["root", "A", "B"],                            vec![("root", "A"), ("root", "B")]),
    ("root", "A"),
    (vec![("root", "B")],                               vec![],                                  vec!["B"]),
  )]
  #[trace]
  fn test_graph_ops_collapse_edge_merges_target_into_source(
    #[case] (nodes, edges): (Vec<&'static str>, Vec<EdgeNames>),
    #[case] collapsed: EdgeNames,
    #[case] (expected_source_outbound, expected_new_edges, expected_children): (Vec<EdgeNames>, Vec<EdgeNames>, Vec<&str>),
  ) -> Result<(), Report> {
    let mut tree = NamedGraph::new(&nodes, &edges)?;
    let (source, target) = collapsed;
    let collapsed_key = tree.edge(source, target)?;
    let expected_source_outbound = edge_keys(&tree, &expected_source_outbound)?;
    let expected_new_edges = edge_keys(&tree, &expected_new_edges)?;

    let (removed_node, removed_edge, new_edges) = tree.graph.collapse_edge(collapsed_key)?;

    let source_outbound = tree.graph.get_node(tree.key(source)).expect("source survives").outbound();
    assert_eq!(expected_source_outbound, source_outbound);
    assert_eq!(expected_new_edges, new_edges);
    assert_eq!(expected_children, tree.children(source));
    assert_eq!(
      (tree.key(target), collapsed_key, false, false),
      (
        removed_node.key(),
        removed_edge.key(),
        tree.graph.get_node(tree.key(target)).is_some(),
        tree.graph.get_edge(collapsed_key).is_some(),
      )
    );
    Ok(())
  }

  #[test]
  fn test_graph_ops_collapse_edge_redirects_other_inbound_edges_to_source() -> Result<(), Report> {
    let mut dag = NamedGraph::new(&["s1", "s2", "t", "leaf"], &[("s1", "t"), ("s2", "t"), ("t", "leaf")])?;
    let collapsed_key = dag.edge("s1", "t")?;
    let other_inbound = dag.edge("s2", "t")?;
    let outbound = dag.edge("t", "leaf")?;

    dag.graph.collapse_edge(collapsed_key)?;

    let s1 = dag.graph.get_node(dag.key("s1")).expect("source survives");
    let s2 = dag.graph.get_node(dag.key("s2")).expect("other parent survives");
    let other_inbound_edge = dag.graph.get_edge(other_inbound).expect("edge survives");
    let outbound_edge = dag.graph.get_edge(outbound).expect("edge survives");
    assert_eq!(
      (
        vec![other_inbound],
        vec![outbound],
        vec![other_inbound],
        ("s2", "s1"),
        ("s1", "leaf"),
      ),
      (
        s1.inbound().to_vec(),
        s1.outbound().to_vec(),
        s2.outbound().to_vec(),
        (
          dag.name(other_inbound_edge.source()),
          dag.name(other_inbound_edge.target())
        ),
        (dag.name(outbound_edge.source()), dag.name(outbound_edge.target())),
      )
    );
    Ok(())
  }

  #[test]
  fn test_graph_ops_collapse_edge_rejects_unknown_edge() -> Result<(), Report> {
    let mut tree = NamedGraph::new(&["root", "A"], &[("root", "A")])?;

    let result = tree.graph.collapse_edge(GraphEdgeKey(9999));

    assert_error!(
      result,
      "Edge 9999 not found. This is an internal error. Please report it to developers."
    );
    Ok(())
  }

  mod helpers {
    use super::super::super::test_graph_support::tests::NamedGraph;
    use crate::edge::GraphEdgeKey;
    use eyre::Report;
    use itertools::Itertools;

    pub(super) type EdgeNames = (&'static str, &'static str);

    pub(super) fn edge_keys(tree: &NamedGraph, edges: &[EdgeNames]) -> Result<Vec<GraphEdgeKey>, Report> {
      edges
        .iter()
        .map(|(source, target)| tree.edge(source, target))
        .try_collect()
    }
  }
}
