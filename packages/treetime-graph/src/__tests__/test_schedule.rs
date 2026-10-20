#[cfg(test)]
mod tests {
  use self::helpers::names;
  use super::super::test_graph_support::tests::helpers::NamedGraph;
  use crate::graph::Graph;
  use crate::schedule::GraphSchedule;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[test]
  fn test_schedule_orders_follow_graph_child_order() -> Result<(), Report> {
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

    let schedule = GraphSchedule::new(&tree.graph)?;

    assert_eq!(
      (
        vec!["root", "CD", "C", "D", "AB", "B", "A"],
        vec!["C", "D", "CD", "B", "A", "AB", "root"],
        vec!["root"],
        vec!["C", "D", "B", "A"],
      ),
      (
        names(&tree, schedule.preorder()),
        names(&tree, schedule.postorder()),
        names(&tree, schedule.roots()),
        names(&tree, schedule.leaves()),
      )
    );
    Ok(())
  }

  #[test]
  fn test_schedule_postorder_lists_a_shared_child_once_before_both_parents() -> Result<(), Report> {
    let dag = NamedGraph::new(&["r", "a", "b"], &[("r", "a"), ("r", "b"), ("a", "b")])?;

    let schedule = GraphSchedule::new(&dag.graph)?;

    assert_eq!(
      (vec!["r", "a", "b"], vec!["b", "a", "r"]),
      (names(&dag, schedule.preorder()), names(&dag, schedule.postorder()))
    );
    Ok(())
  }

  #[test]
  fn test_schedule_preorder_runs_a_shared_child_after_its_last_parent() -> Result<(), Report> {
    let dag = NamedGraph::new(
      &["r", "a", "b", "c", "s"],
      &[("r", "a"), ("r", "b"), ("a", "s"), ("b", "c"), ("c", "s")],
    )?;

    let schedule = GraphSchedule::new(&dag.graph)?;

    assert_eq!(
      (vec!["r", "a", "b", "c", "s"], vec!["s", "a", "c", "b", "r"]),
      (names(&dag, schedule.preorder()), names(&dag, schedule.postorder()))
    );
    Ok(())
  }

  #[test]
  fn test_schedule_slots_give_the_position_of_an_edge_at_each_end() -> Result<(), Report> {
    let dag = NamedGraph::new(&["p", "q", "x", "s"], &[("p", "x"), ("q", "s"), ("p", "s"), ("x", "s")])?;

    let schedule = GraphSchedule::new(&dag.graph)?;

    let (p_s, q_s, x_s) = (dag.edge("p", "s")?, dag.edge("q", "s")?, dag.edge("x", "s")?);
    assert_eq!(
      (
        vec![("q", q_s), ("p", p_s), ("x", x_s)],
        vec![("x", dag.edge("p", "x")?), ("s", p_s)],
        (1, 0, 0),
        (1, 0, 2),
        vec!["p", "q"],
        vec!["s"],
      ),
      (
        schedule
          .parents(dag.key("s"))
          .iter()
          .map(|&(parent, edge)| (dag.name(parent), edge))
          .collect_vec(),
        schedule
          .children(dag.key("p"))
          .iter()
          .map(|&(child, edge)| (dag.name(child), edge))
          .collect_vec(),
        (
          schedule.child_slot(p_s),
          schedule.child_slot(x_s),
          schedule.child_slot(q_s)
        ),
        (
          schedule.parent_slot(p_s),
          schedule.parent_slot(q_s),
          schedule.parent_slot(x_s)
        ),
        names(&dag, schedule.roots()),
        names(&dag, schedule.leaves()),
      )
    );
    Ok(())
  }

  #[test]
  fn test_schedule_single_node_is_root_and_leaf() -> Result<(), Report> {
    let tree = NamedGraph::new(&["root"], &[])?;

    let schedule = GraphSchedule::new(&tree.graph)?;

    assert_eq!(
      (vec!["root"], vec!["root"], vec!["root"], vec!["root"]),
      (
        names(&tree, schedule.preorder()),
        names(&tree, schedule.postorder()),
        names(&tree, schedule.roots()),
        names(&tree, schedule.leaves()),
      )
    );
    Ok(())
  }

  #[test]
  fn test_schedule_skips_deleted_key_gaps() -> Result<(), Report> {
    let mut graph = Graph::new();
    let root = graph.add_node();
    let removed = graph.add_node();
    let child = graph.add_node();
    graph.add_edge(root, child)?;
    graph.remove_node(removed)?;
    graph.build()?;

    let schedule = GraphSchedule::new(&graph)?;

    assert_eq!(
      (vec![root, child], vec![child, root], 0, 0),
      (
        schedule.preorder().to_vec(),
        schedule.postorder().to_vec(),
        schedule.parents(removed).len(),
        schedule.children(removed).len(),
      )
    );
    Ok(())
  }

  #[test]
  fn test_schedule_rejects_graph_without_root() -> Result<(), Report> {
    let cycle = NamedGraph::new(&["A", "B"], &[("A", "B"), ("B", "A")])?;

    assert_error!(
      GraphSchedule::new(&cycle.graph).map(|_| ()),
      "The graph has no root: every node has a parent"
    );
    Ok(())
  }

  #[test]
  fn test_schedule_rejects_empty_graph() -> Result<(), Report> {
    let graph = Graph::new();

    assert_error!(
      GraphSchedule::new(&graph).map(|_| ()),
      "The graph has no root: every node has a parent"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::reached_by_the_root(   &["r", "a", "b"],                &[("r", "a"), ("a", "b"), ("b", "a")],             "a")]
  #[case::not_reached(           &["root", "A", "X", "Y"],        &[("root", "A"), ("X", "Y"), ("Y", "X")],          "X")]
  #[case::node_below_the_cycle(  &["root", "A", "Z", "X", "Y"],   &[("root", "A"), ("X", "Y"), ("Y", "X"), ("Y", "Z")], "Y")]
  #[trace]
  fn test_schedule_rejects_cycle(
    #[case] nodes: &[&'static str],
    #[case] edges: &[(&'static str, &'static str)],
    #[case] on_cycle: &str,
  ) -> Result<(), Report> {
    let graph = NamedGraph::new(nodes, edges)?;

    assert_error!(
      GraphSchedule::new(&graph.graph).map(|_| ()),
      format!("The graph has a cycle through node {}", graph.key(on_cycle))
    );
    Ok(())
  }

  #[test]
  fn test_schedule_rejects_node_that_is_its_own_parent() -> Result<(), Report> {
    let mut graph = Graph::new();
    let root = graph.add_node();
    let node = graph.add_node();
    let other = graph.add_node();
    graph.add_edge(root, node)?;
    let edge = graph.add_edge(node, other)?;
    graph.edge_mut(edge).set_target(node);
    graph.node_mut(other).inbound_mut().clear();
    graph.node_mut(node).inbound_mut().push(edge);
    graph.build()?;

    assert_error!(
      GraphSchedule::new(&graph).map(|_| ()),
      format!("The graph has a cycle through node {node}")
    );
    Ok(())
  }

  #[test]
  fn test_schedule_require_tree_rejects_two_roots() -> Result<(), Report> {
    let forest = NamedGraph::new(&["r1", "r2", "A", "B"], &[("r1", "A"), ("r2", "B")])?;
    let schedule = GraphSchedule::new(&forest.graph)?;

    let expected = format!(
      "The graph is not a tree: it has 2 roots (node {}, node {}), but a tree has one",
      forest.key("r1"),
      forest.key("r2")
    );
    assert_error!(schedule.require_tree(), expected);
    assert_error!(schedule.into_tree().map(|_| ()), expected);
    Ok(())
  }

  #[test]
  fn test_schedule_require_tree_rejects_node_with_two_parents() -> Result<(), Report> {
    let network = NamedGraph::new(
      &["root", "A", "B", "C"],
      &[("root", "A"), ("root", "B"), ("A", "C"), ("B", "C")],
    )?;
    let schedule = GraphSchedule::new(&network.graph)?;

    let expected = format!(
      "The graph is not a tree: node {} has 2 parents, but a node of a tree has at most one",
      network.key("C")
    );
    assert_error!(schedule.require_tree(), expected);
    assert_error!(schedule.into_tree().map(|_| ()), expected);
    Ok(())
  }

  #[test]
  fn test_schedule_tree_gives_root_parent_and_children() -> Result<(), Report> {
    let tree = NamedGraph::new(
      &["root", "AB", "A", "B", "C"],
      &[("root", "AB"), ("root", "C"), ("AB", "A"), ("AB", "B")],
    )?;

    let tree_schedule = GraphSchedule::new(&tree.graph)?.into_tree()?;

    assert_eq!(
      (
        "root",
        None,
        Some((tree.key("AB"), tree.edge("AB", "B")?)),
        vec![("AB", tree.edge("root", "AB")?), ("C", tree.edge("root", "C")?)],
        0,
        vec!["root", "AB", "A", "B", "C"],
        vec!["A", "B", "AB", "C", "root"],
      ),
      (
        tree.name(tree_schedule.root()),
        tree_schedule.parent(tree.key("root")),
        tree_schedule.parent(tree.key("B")),
        tree_schedule
          .children(tree.key("root"))
          .iter()
          .map(|&(child, edge)| (tree.name(child), edge))
          .collect_vec(),
        tree_schedule.children(tree.key("A")).len(),
        names(&tree, tree_schedule.preorder()),
        names(&tree, tree_schedule.schedule().postorder()),
      )
    );
    Ok(())
  }

  mod helpers {
    use super::super::super::test_graph_support::tests::helpers::NamedGraph;
    use crate::node::GraphNodeKey;
    use itertools::Itertools;

    pub(super) fn names(graph: &NamedGraph, keys: &[GraphNodeKey]) -> Vec<&'static str> {
      keys.iter().map(|&key| graph.name(key)).collect_vec()
    }
  }
}
