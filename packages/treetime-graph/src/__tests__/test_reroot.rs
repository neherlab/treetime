#[cfg(test)]
mod tests {
  use super::super::test_graph_support::tests::NamedGraph;
  use crate::reroot::remove_stem_root;
  use eyre::Report;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_remove_stem_root_makes_the_only_child_the_root() -> Result<(), Report> {
    let mut tree = NamedGraph::new(
      &["stem", "root", "A", "B"],
      &[("stem", "root"), ("root", "A"), ("root", "B")],
    )?;
    let stem_edge = tree.edge("stem", "root")?;

    let stem_key = tree.key("stem");
    let stem = remove_stem_root(&mut tree.graph, stem_key)?.expect("a one-child root is a stem");

    assert_eq!(
      (tree.key("stem"), stem_edge, tree.key("root")),
      (stem.removed_node_key, stem.removed_edge_key, stem.new_root_key)
    );
    assert!(tree.graph.get_node(tree.key("stem")).is_none());
    assert!(tree.graph.get_edge(stem_edge).is_none());
    assert_eq!(tree.key("root"), tree.graph.get_exactly_one_root()?.key());
    assert_eq!(vec!["A", "B"], tree.children("root"));
    Ok(())
  }

  #[test]
  fn test_remove_stem_root_keeps_a_root_with_several_children() -> Result<(), Report> {
    let mut tree = NamedGraph::new(&["root", "A", "B"], &[("root", "A"), ("root", "B")])?;

    let root_key = tree.key("root");
    let stem = remove_stem_root(&mut tree.graph, root_key)?;

    assert!(stem.is_none());
    assert_eq!(3, tree.graph.get_nodes().count());
    assert_eq!(tree.key("root"), tree.graph.get_exactly_one_root()?.key());
    Ok(())
  }
}
