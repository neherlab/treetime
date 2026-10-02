#[cfg(test)]
mod tests {
  use crate::assign_node_names::restrict_node_names;
  use crate::graph::Graph;
  use crate::node::GraphNodeKey;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_utils::o;

  #[test]
  fn test_restrict_node_names_keeps_graph_nodes_only() -> Result<(), Report> {
    let mut graph = Graph::new();
    let root = graph.add_node();
    let named = graph.add_node();
    let unnamed = graph.add_node();
    let without_entry = graph.add_node();
    graph.add_edge(root, named)?;
    graph.add_edge(root, unnamed)?;
    graph.add_edge(root, without_entry)?;
    graph.build()?;
    let removed = GraphNodeKey(99);

    let names = btreemap! {
      root => Some(o!("root")),
      named => Some(o!("A")),
      unnamed => None,
      removed => Some(o!("gone")),
    };

    let actual = restrict_node_names(&names, &graph);

    assert_eq!(
      actual,
      btreemap! {
        root => Some(o!("root")),
        named => Some(o!("A")),
        unnamed => None,
        without_entry => None,
      }
    );
    Ok(())
  }
}
