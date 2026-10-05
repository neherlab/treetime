#[cfg(test)]
mod tests {
  use self::helpers::{graph_from_parents, names};
  use super::super::test_graph_support::tests::NamedGraph;
  use crate::tree_view::TreeView;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use treetime_utils::assert_error;

  #[test]
  fn test_tree_view_preorder_lists_parents_before_children_in_graph_child_order() -> Result<(), Report> {
    let tree = NamedGraph::new(
      &["root", "AB", "CD", "A", "B", "C", "D"],
      &[
        ("root", "CD"),
        ("root", "AB"),
        ("AB", "B"),
        ("AB", "A"),
        ("CD", "C"),
        ("CD", "D"),
      ],
    )?;

    let view = TreeView::new(&tree.graph)?;

    assert_eq!(tree.key("root"), view.root());
    assert_eq!(
      vec!["root", "CD", "C", "D", "AB", "B", "A"],
      names(&tree, view.preorder())
    );
    Ok(())
  }

  #[test]
  fn test_tree_view_parent_and_children_follow_edges() -> Result<(), Report> {
    let tree = NamedGraph::new(
      &["root", "AB", "A", "B", "C"],
      &[("root", "AB"), ("root", "C"), ("AB", "A"), ("AB", "B")],
    )?;

    let view = TreeView::new(&tree.graph)?;

    assert_eq!(None, view.parent(tree.key("root")));
    assert_eq!(
      Some((tree.key("AB"), tree.edge("AB", "B")?)),
      view.parent(tree.key("B"))
    );
    let children = view
      .children(tree.key("root"))
      .iter()
      .map(|&(child, edge)| (tree.name(child), edge))
      .collect_vec();
    assert_eq!(
      vec![("AB", tree.edge("root", "AB")?), ("C", tree.edge("root", "C")?)],
      children
    );
    assert_eq!(0, view.children(tree.key("A")).len());
    Ok(())
  }

  #[test]
  fn test_tree_view_single_node_is_root_without_children() -> Result<(), Report> {
    let tree = NamedGraph::new(&["root"], &[])?;

    let view = TreeView::new(&tree.graph)?;

    assert_eq!(vec!["root"], names(&tree, view.preorder()));
    assert_eq!(None, view.parent(tree.key("root")));
    assert_eq!(0, view.children(tree.key("root")).len());
    Ok(())
  }

  #[test]
  fn test_tree_view_rejects_two_roots() -> Result<(), Report> {
    let forest = NamedGraph::new(&["r1", "r2", "A", "B"], &[("r1", "A"), ("r2", "B")])?;

    assert_error!(
      TreeView::new(&forest.graph),
      format!(
        "The graph is not a tree: it has 2 roots (node {}, node {}), but a tree has one",
        forest.key("r1"),
        forest.key("r2")
      )
    );
    Ok(())
  }

  #[test]
  fn test_tree_view_rejects_graph_without_root() -> Result<(), Report> {
    let cycle = NamedGraph::new(&["A", "B"], &[("A", "B"), ("B", "A")])?;

    assert_error!(
      TreeView::new(&cycle.graph),
      "The graph is not a tree: every node has a parent, so there is no root"
    );
    Ok(())
  }

  #[test]
  fn test_tree_view_rejects_node_with_two_parents() -> Result<(), Report> {
    let network = NamedGraph::new(
      &["root", "A", "B", "C"],
      &[("root", "A"), ("root", "B"), ("A", "C"), ("B", "C")],
    )?;

    assert_error!(
      TreeView::new(&network.graph),
      format!(
        "The graph is not a tree: node {} has 2 parents, but a node of a tree has at most one",
        network.key("C")
      )
    );
    Ok(())
  }

  #[test]
  fn test_tree_view_rejects_cycle_that_the_root_does_not_reach() -> Result<(), Report> {
    let graph = NamedGraph::new(&["root", "A", "X", "Y"], &[("root", "A"), ("X", "Y"), ("Y", "X")])?;

    assert_error!(
      TreeView::new(&graph.graph),
      format!(
        "The graph is not a tree: node {} is on a cycle, which the root does not reach",
        graph.key("X")
      )
    );
    Ok(())
  }

  #[test]
  fn test_tree_view_rejects_unreachable_node_below_a_cycle() -> Result<(), Report> {
    let graph = NamedGraph::new(
      &["root", "A", "Z", "X", "Y"],
      &[("root", "A"), ("X", "Y"), ("Y", "X"), ("Y", "Z")],
    )?;

    assert_error!(
      TreeView::new(&graph.graph),
      format!(
        "The graph is not a tree: node {} is on a cycle, which the root does not reach",
        graph.key("Y")
      )
    );
    Ok(())
  }

  proptest! {
    #[test]
    fn test_prop_tree_view_matches_graph_traversal_and_edges(
      parents in (1_usize..40).prop_flat_map(|count| {
        (0..count).map(|index| 0..index.max(1)).collect_vec()
      }),
    ) {
      let graph = graph_from_parents(&parents).unwrap();

      let view = TreeView::new(&graph).unwrap();

      let mut expected = vec![];
      graph.iter_depth_first_preorder_forward(|node| {
        expected.push(node.key);
        Ok(())
      }).unwrap();
      prop_assert_eq!(&expected, &view.preorder().to_vec());
      for node in graph.get_nodes() {
        let key = node.key();
        prop_assert_eq!(graph.node_parent(key).unwrap(), view.parent(key));
        prop_assert_eq!(graph.children_keys_of(node).collect_vec(), view.children(key).to_vec());
        if let Some((parent, _)) = view.parent(key) {
          let position = |target| view.preorder().iter().position(|&key| key == target);
          prop_assert!(position(parent) < position(key));
        }
      }
    }
  }

  mod helpers {
    use super::super::super::test_graph_support::tests::NamedGraph;
    use crate::graph::Graph;
    use crate::node::GraphNodeKey;
    use eyre::Report;
    use itertools::Itertools;

    pub(super) fn names(tree: &NamedGraph, keys: &[GraphNodeKey]) -> Vec<&'static str> {
      keys.iter().map(|&key| tree.name(key)).collect_vec()
    }

    pub(super) fn graph_from_parents(parents: &[usize]) -> Result<Graph, Report> {
      let mut graph = Graph::new();
      let keys = std::iter::repeat_with(|| graph.add_node())
        .take(parents.len())
        .collect_vec();
      for (child, &parent) in parents.iter().enumerate().skip(1) {
        graph.add_edge(keys[parent], keys[child])?;
      }
      graph.build()?;
      Ok(graph)
    }
  }
}
