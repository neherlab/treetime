#[cfg(test)]
mod tests {
  use self::helpers::{matches_by_name, named_tree, tips_by_name};
  use crate::clades::{CladeNode, TipIds, clade_keys, first_node_by_key, graph_clade_keys};
  use eyre::Report;
  use itertools::Itertools;
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_utils::assert_error;

  const NODES: &[&str] = &["root", "AB", "A", "B", "C"];
  const EDGES: &[(&str, &str)] = &[("root", "AB"), ("root", "C"), ("AB", "A"), ("AB", "B")];

  #[test]
  fn test_clades_graph_keys_hold_the_tips_below_each_node() -> Result<(), Report> {
    let (graph, names) = named_tree(NODES, EDGES)?;
    let mut ids = TipIds::default();

    let keys = graph_clade_keys(&graph, &names, &mut ids)?;

    let expected = btreemap! {
      "root" => btreeset! {"A", "B", "C"},
      "AB"   => btreeset! {"A", "B"},
      "A"    => btreeset! {"A"},
      "B"    => btreeset! {"B"},
      "C"    => btreeset! {"C"},
    };
    assert_eq!(expected, tips_by_name(&keys, &names, &ids));
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::reordered_children(
    &["C", "B", "root", "A", "AB"],
    &[("root", "C"), ("AB", "B"), ("root", "AB"), ("AB", "A")],
    btreemap! {"root" => Some("root"), "AB" => Some("AB"), "A" => Some("A"), "B" => Some("B"), "C" => Some("C")},
  )]
  #[case::different_clades(
    &["root", "AC", "A", "B", "C"],
    &[("root", "AC"), ("root", "B"), ("AC", "A"), ("AC", "C")],
    btreemap! {"root" => Some("root"), "AB" => None,       "A" => Some("A"), "B" => Some("B"), "C" => Some("C")},
  )]
  #[trace]
  fn test_clades_keys_of_two_trees_compare_through_one_interner(
    #[case] other_nodes: &[&str],
    #[case] other_edges: &[(&str, &str)],
    #[case] expected: BTreeMap<&str, Option<&str>>,
  ) -> Result<(), Report> {
    let (graph, names) = named_tree(NODES, EDGES)?;
    let (other_graph, other_names) = named_tree(other_nodes, other_edges)?;
    let mut ids = TipIds::default();

    let keys = graph_clade_keys(&graph, &names, &mut ids)?;
    let other_keys = graph_clade_keys(&other_graph, &other_names, &mut ids)?;

    assert_eq!(expected, matches_by_name(&keys, &names, &other_keys, &other_names));
    Ok(())
  }

  #[test]
  fn test_clades_first_node_by_key_takes_the_first_node_in_preorder() {
    let nodes = [
      CladeNode {
        parent: None,
        tip: None,
      },
      CladeNode {
        parent: Some(0),
        tip: None,
      },
      CladeNode {
        parent: Some(1),
        tip: Some("A"),
      },
      CladeNode {
        parent: Some(0),
        tip: Some("B"),
      },
    ];
    let mut ids = TipIds::default();

    let keys = clade_keys(nodes, &mut ids);
    let first = first_node_by_key(keys.iter().map(Vec::as_slice).enumerate());

    let expected_keys = vec![vec![0, 1], vec![0], vec![0], vec![1]];
    let expected_first = btreemap! {
      [0, 1].as_slice() => 0,
      [0].as_slice()    => 1,
      [1].as_slice()    => 3,
    };
    assert_eq!((&expected_keys, expected_first), (&keys, first));
  }

  #[test]
  fn test_clades_keys_keep_a_tip_once_per_occurrence() {
    let nodes = [
      CladeNode {
        parent: None,
        tip: None,
      },
      CladeNode {
        parent: Some(0),
        tip: Some("A"),
      },
      CladeNode {
        parent: Some(0),
        tip: Some("A"),
      },
    ];
    let mut ids = TipIds::default();

    let keys = clade_keys(nodes, &mut ids);

    assert_eq!(vec![vec![0, 0], vec![0], vec![0]], keys);
  }

  #[test]
  #[should_panic(expected = "node 1 has parent 2, which does not precede it in preorder")]
  fn test_clades_keys_reject_nodes_out_of_preorder() {
    let nodes = [
      CladeNode {
        parent: None,
        tip: None,
      },
      CladeNode {
        parent: Some(2),
        tip: Some("A"),
      },
      CladeNode {
        parent: Some(0),
        tip: None,
      },
    ];
    clade_keys(nodes, &mut TipIds::default());
  }

  #[test]
  fn test_clades_graph_keys_reject_a_node_with_two_parents() -> Result<(), Report> {
    let (graph, names) = named_tree(
      &["root", "X", "Y", "A"],
      &[("root", "X"), ("root", "Y"), ("X", "A"), ("Y", "A")],
    )?;

    assert_error!(
      graph_clade_keys(&graph, &names, &mut TipIds::default()),
      "The graph is not a tree: node 3 has 2 parents, but a node of a tree has at most one"
    );
    Ok(())
  }

  #[test]
  fn test_clades_tip_ids_give_back_the_interned_names() {
    let mut ids = TipIds::default();

    let interned = ["B", "A", "B"].map(|name| ids.id(name));

    assert_eq!(
      ([0, 1, 0], vec!["B", "A", "B"]),
      (interned, interned.iter().map(|&id| ids.name(id)).collect_vec())
    );
  }

  mod helpers {
    use crate::clades::{TipIds, first_node_by_key};
    use crate::graph::Graph;
    use crate::node::GraphNodeKey;
    use crate::node_names::NodeNames;
    use eyre::Report;
    use indexmap::IndexMap;
    use std::collections::{BTreeMap, BTreeSet};

    pub(super) fn named_tree(nodes: &[&str], edges: &[(&str, &str)]) -> Result<(Graph, NodeNames), Report> {
      let mut graph = Graph::new();
      let keys: BTreeMap<_, _> = nodes.iter().map(|&name| (name, graph.add_node())).collect();
      for (source, target) in edges {
        graph.add_edge(keys[source], keys[target])?;
      }
      graph.build()?;
      let input = keys.iter().map(|(&name, &key)| (key, Some(name.to_owned()))).collect();
      let names = NodeNames::from_input(input, &graph)?.names;
      Ok((graph, names))
    }

    pub(super) fn tips_by_name<'a>(
      keys: &IndexMap<GraphNodeKey, Vec<usize>>,
      names: &'a NodeNames,
      ids: &TipIds<'a>,
    ) -> BTreeMap<&'a str, BTreeSet<&'a str>> {
      keys
        .iter()
        .map(|(&node, key)| (&names[node], key.iter().map(|&id| ids.name(id)).collect()))
        .collect()
    }

    pub(super) fn matches_by_name<'a>(
      keys: &IndexMap<GraphNodeKey, Vec<usize>>,
      names: &'a NodeNames,
      other_keys: &IndexMap<GraphNodeKey, Vec<usize>>,
      other_names: &'a NodeNames,
    ) -> BTreeMap<&'a str, Option<&'a str>> {
      let other_nodes = first_node_by_key(other_keys.iter().map(|(&node, key)| (node, key.as_slice())));
      keys
        .iter()
        .map(|(&node, key)| {
          let other = other_nodes.get(key.as_slice()).map(|&other| &other_names[other]);
          (&names[node], other)
        })
        .collect()
    }
  }
}
