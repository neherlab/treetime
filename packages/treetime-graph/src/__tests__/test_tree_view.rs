#[cfg(test)]
mod tests {
  use self::helpers::names;
  use super::super::test_graph_support::tests::helpers::{NamedGraph, graph_from_parents};
  use crate::schedule::GraphSchedule;
  use crate::tree_view::TreeView;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;

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

  proptest! {
    #[test]
    fn test_prop_tree_view_matches_graph_traversal_and_edges(
      parents in (1_usize..40).prop_flat_map(|count| {
        (0..count).map(|index| 0..index.max(1)).collect_vec()
      }),
    ) {
      let graph = graph_from_parents(&parents).unwrap();

      let view = TreeView::new(&graph).unwrap();

      let expected = GraphSchedule::new(&graph).unwrap().preorder().to_vec();
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
    use super::super::super::test_graph_support::tests::helpers::NamedGraph;
    use crate::node::GraphNodeKey;
    use itertools::Itertools;

    pub(super) fn names(tree: &NamedGraph, keys: &[GraphNodeKey]) -> Vec<&'static str> {
      keys.iter().map(|&key| tree.name(key)).collect_vec()
    }
  }
}
