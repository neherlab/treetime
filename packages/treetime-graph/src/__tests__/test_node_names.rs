#[cfg(test)]
mod tests {
  use self::helpers::{names_by_key, tree};
  use crate::node::GraphNodeKey;
  use crate::node_names::{InputNodeNames, NodeNames};
  use eyre::Report;
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_utils::{o, vec_of_owned};

  #[test]
  fn test_node_names_from_input_names_unnamed_nodes_in_preorder() -> Result<(), Report> {
    let graph = tree()?;
    let input = btreemap! {
      GraphNodeKey(2) => Some(o!("A")),
      GraphNodeKey(3) => Some(o!("B")),
      GraphNodeKey(4) => Some(o!("X")),
      GraphNodeKey(5) => Some(o!("C")),
    };

    let InputNodeNames {
      names,
      unnamed_nodes,
      duplicate_names,
    } = NodeNames::from_input(input, &graph)?;

    assert_eq!(
      vec_of_owned!["NODE_0000000", "NODE_0000001", "A", "B", "X", "C"],
      names_by_key(&graph, &names)
    );
    assert_eq!(btreeset! { GraphNodeKey(0), GraphNodeKey(1) }, unnamed_nodes);
    assert_eq!(Vec::<String>::new(), duplicate_names);
    Ok(())
  }

  #[test]
  fn test_node_names_from_input_skips_generated_names_taken_by_the_input() -> Result<(), Report> {
    let graph = tree()?;
    let input = btreemap! {
      GraphNodeKey(1) => Some(o!("")),
      GraphNodeKey(2) => Some(o!("NODE_0000000")),
      GraphNodeKey(3) => Some(o!("B")),
      GraphNodeKey(4) => None,
      GraphNodeKey(5) => Some(o!("NODE_0000002")),
    };

    let actual = NodeNames::from_input(input, &graph)?.names;

    assert_eq!(
      vec_of_owned![
        "NODE_0000001",
        "NODE_0000003",
        "NODE_0000000",
        "B",
        "NODE_0000004",
        "NODE_0000002"
      ],
      names_by_key(&graph, &actual)
    );
    Ok(())
  }

  #[test]
  fn test_node_names_from_input_reports_repeated_input_names() -> Result<(), Report> {
    let graph = tree()?;
    let input: BTreeMap<_, _> = (0..6).map(|key| (GraphNodeKey(key), Some(o!("A")))).collect();

    let actual = NodeNames::from_input(input, &graph)?;

    assert_eq!(vec_of_owned!["A"], actual.duplicate_names);
    assert_eq!(
      vec_of_owned!["A", "A", "A", "A", "A", "A"],
      names_by_key(&graph, &actual.names)
    );
    Ok(())
  }

  #[test]
  fn test_node_names_name_new_node_continues_after_every_name_in_use() -> Result<(), Report> {
    let graph = tree()?;
    let input = btreemap! {
      GraphNodeKey(2) => Some(o!("A")),
      GraphNodeKey(3) => Some(o!("NODE_0000002")),
      GraphNodeKey(4) => Some(o!("NODE_0000004")),
      GraphNodeKey(5) => Some(o!("C")),
    };
    let mut names = NodeNames::from_input(input, &graph)?.names;

    names.name_new_node(GraphNodeKey(6));
    names.name_new_node(GraphNodeKey(7));
    names.name_new_node(GraphNodeKey(8));

    let actual = [6, 7, 8].map(|key| names[GraphNodeKey(key)].to_owned());
    assert_eq!(["NODE_0000003", "NODE_0000005", "NODE_0000006"], actual);
    Ok(())
  }

  mod helpers {
    use crate::__tests__::test_graph_support::tests::helpers::graph_from_parents;
    use crate::graph::Graph;
    use crate::node_names::NodeNames;
    use eyre::Report;

    pub(super) fn tree() -> Result<Graph, Report> {
      graph_from_parents(&[0, 0, 1, 1, 0, 4])
    }

    pub(super) fn names_by_key(graph: &Graph, names: &NodeNames) -> Vec<String> {
      graph.get_nodes().map(|node| names[node.key()].to_owned()).collect()
    }
  }
}
