#[cfg(test)]
mod tests {
  use self::helpers::{LengthEdge, Tree};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::leaf_of_multifurcation(
    ("R", vec![("R", "A", Some(1.0)), ("R", "B", Some(2.0)), ("R", "C", Some(3.0))]),
    vec!["A"],
    ("R", vec![("R", "B", Some(2.0)), ("R", "C", Some(3.0))]),
  )]
  #[case::cherry_parent_merges_into_child(
    ("R", vec![("R", "X", Some(3.0)), ("R", "C", Some(4.0)), ("X", "A", Some(1.0)), ("X", "B", Some(2.0))]),
    vec!["A"],
    ("R", vec![("R", "B", Some(5.0)), ("R", "C", Some(4.0))]),
  )]
  #[case::merge_keeps_known_length(
    ("R", vec![("R", "X", None), ("R", "C", Some(4.0)), ("X", "A", Some(1.0)), ("X", "B", Some(2.0))]),
    vec!["A"],
    ("R", vec![("R", "B", Some(2.0)), ("R", "C", Some(4.0))]),
  )]
  #[case::whole_cherry(
    ("R", vec![("R", "X", Some(3.0)), ("R", "C", Some(4.0)), ("R", "D", Some(5.0)), ("X", "A", Some(1.0)), ("X", "B", Some(2.0))]),
    vec!["A", "B"],
    ("R", vec![("R", "C", Some(4.0)), ("R", "D", Some(5.0))]),
  )]
  #[case::root_with_one_child_left(
    ("R", vec![("R", "A", Some(1.0)), ("R", "X", Some(2.0)), ("X", "B", Some(1.0)), ("X", "C", Some(1.5))]),
    vec!["A"],
    ("X", vec![("X", "B", Some(1.0)), ("X", "C", Some(1.5))]),
  )]
  #[case::nested_parents(
    ("R", vec![("R", "X", Some(1.0)), ("R", "D", Some(9.0)), ("X", "Y", Some(2.0)), ("X", "C", Some(3.0)), ("Y", "A", Some(4.0)), ("Y", "B", Some(5.0))]),
    vec!["A", "C"],
    ("R", vec![("R", "B", Some(8.0)), ("R", "D", Some(9.0))]),
  )]
  #[trace]
  fn test_remove_leaves_merges_parents_left_with_one_child(
    #[case] (root, edges): (&'static str, Vec<LengthEdge>),
    #[case] removed: Vec<&'static str>,
    #[case] expected: (&'static str, Vec<LengthEdge>),
  ) -> Result<(), Report> {
    let mut tree = Tree::new(root, &edges)?;

    tree.remove(&removed)?;

    assert_eq!(expected, (tree.root(), tree.edges()));
    Ok(())
  }

  #[test]
  fn test_remove_leaves_rejects_removing_every_leaf() -> Result<(), Report> {
    let mut tree = Tree::new("R", &[("R", "A", Some(1.0)), ("R", "B", Some(2.0))])?;

    assert_error!(tree.remove(&["A", "B"]), "Removing the tips would leave an empty tree");
    Ok(())
  }

  #[test]
  fn test_remove_leaves_rejects_internal_node() -> Result<(), Report> {
    let mut tree = Tree::new(
      "R",
      &[
        ("R", "X", Some(1.0)),
        ("R", "C", Some(2.0)),
        ("X", "A", Some(1.0)),
        ("X", "B", Some(1.0)),
      ],
    )?;

    assert_error!(
      tree.remove(&["X"]),
      "Node 1 to remove is not a leaf. This is an internal error. Please report it to developers."
    );
    Ok(())
  }

  mod helpers {
    use crate::edge::GraphEdgeKey;
    use crate::graph::Graph;
    use crate::key_map::EdgeMap;
    use crate::node::GraphNodeKey;
    use crate::remove_leaves::remove_leaves;
    use eyre::Report;
    use itertools::Itertools;
    use std::collections::{BTreeMap, BTreeSet};
    use std::iter;

    pub(crate) type LengthEdge = (&'static str, &'static str, Option<f64>);

    pub(crate) struct Tree {
      graph: Graph,
      branch_lengths: EdgeMap<Option<f64>>,
      keys: BTreeMap<&'static str, GraphNodeKey>,
    }

    impl Tree {
      pub(crate) fn new(root: &'static str, edges: &[LengthEdge]) -> Result<Self, Report> {
        let mut graph = Graph::new();
        let names = edges
          .iter()
          .flat_map(|(source, target, _)| [*source, *target])
          .filter(|name| *name != root);
        let keys: BTreeMap<_, _> = iter::once(root)
          .chain(names)
          .unique()
          .map(|name| (name, graph.add_node()))
          .collect();
        let branch_lengths = edges
          .iter()
          .map(
            |(source, target, length)| -> Result<(GraphEdgeKey, Option<f64>), Report> {
              Ok((graph.add_edge(keys[source], keys[target])?, *length))
            },
          )
          .collect::<Result<EdgeMap<_>, Report>>()?;
        graph.build()?;
        Ok(Self {
          graph,
          branch_lengths,
          keys,
        })
      }

      pub(crate) fn remove(&mut self, names: &[&str]) -> Result<(), Report> {
        let leaves: BTreeSet<_> = names.iter().map(|name| self.keys[name]).collect();
        remove_leaves(&mut self.graph, &mut self.branch_lengths, &leaves)
      }

      pub(crate) fn root(&self) -> &'static str {
        self.name(self.graph.root_key().expect("one root"))
      }

      pub(crate) fn edges(&self) -> Vec<LengthEdge> {
        self
          .graph
          .get_edges()
          .map(|edge| {
            (
              self.name(edge.source()),
              self.name(edge.target()),
              self.branch_lengths[&edge.key()],
            )
          })
          .sorted_by_key(|(source, target, _)| (*source, *target))
          .collect()
      }

      fn name(&self, key: GraphNodeKey) -> &'static str {
        self
          .keys
          .iter()
          .find_map(|(name, node_key)| (*node_key == key).then_some(*name))
          .expect("every node has a name")
      }
    }
  }
}
