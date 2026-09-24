#[cfg(test)]
mod tests {
  use super::super::test_graph_support::tests::NamedGraph;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::{assert_error, make_error};

  use self::helpers::{Traversal, fixture_balanced, traverse};

  #[rustfmt::skip]
  #[rstest]
  #[case::depth_first_preorder_forward( Traversal::DepthFirstPreorderForward,  vec!["root", "AB", "A", "B", "CD", "C", "D"])]
  #[case::depth_first_postorder_forward(Traversal::DepthFirstPostorderForward, vec!["A", "B", "AB", "C", "D", "CD", "root"])]
  #[case::breadth_first_forward(        Traversal::BreadthFirstForward,        vec!["root", "AB", "CD", "A", "B", "C", "D"])]
  #[case::breadth_first_backward(       Traversal::BreadthFirstBackward,       vec!["D", "C", "B", "A", "CD", "AB", "root"])]
  #[trace]
  fn test_graph_traverse_visits_every_node_in_order(
    #[case] traversal: Traversal,
    #[case] expected: Vec<&str>,
  ) -> Result<(), Report> {
    let tree = fixture_balanced()?;

    let mut actual = vec![];
    traverse(&tree, traversal, |name| {
      actual.push(name);
      Ok(())
    })?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::depth_first_preorder_forward_at_first_visit( (Traversal::DepthFirstPreorderForward,  "root"), vec!["root"])]
  #[case::depth_first_preorder_forward_at_internal(    (Traversal::DepthFirstPreorderForward,  "AB"),   vec!["root", "AB"])]
  #[case::depth_first_postorder_forward_at_internal(   (Traversal::DepthFirstPostorderForward, "AB"),   vec!["A", "B", "AB"])]
  #[case::depth_first_postorder_forward_at_last_visit( (Traversal::DepthFirstPostorderForward, "root"), vec!["A", "B", "AB", "C", "D", "CD", "root"])]
  #[case::breadth_first_forward_at_first_visit(        (Traversal::BreadthFirstForward,        "root"), vec!["root"])]
  #[case::breadth_first_forward_at_internal(           (Traversal::BreadthFirstForward,        "AB"),   vec!["root", "AB"])]
  #[case::breadth_first_backward_at_internal(          (Traversal::BreadthFirstBackward,       "AB"),   vec!["D", "C", "B", "A", "CD", "AB"])]
  #[case::breadth_first_backward_at_last_visit(        (Traversal::BreadthFirstBackward,       "root"), vec!["D", "C", "B", "A", "CD", "AB", "root"])]
  #[trace]
  fn test_graph_traverse_stops_at_first_error(
    #[case] (traversal, failing): (Traversal, &str),
    #[case] expected_visited: Vec<&str>,
  ) -> Result<(), Report> {
    let tree = fixture_balanced()?;

    let mut visited = vec![];
    let result = traverse(&tree, traversal, |name| {
      visited.push(name);
      if name == failing {
        return make_error!("boom {name}");
      }
      Ok(())
    });

    assert_error!(result, format!("boom {failing}"));
    assert_eq!(expected_visited, visited);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::depth_first_preorder_forward( Traversal::DepthFirstPreorderForward)]
  #[case::depth_first_postorder_forward(Traversal::DepthFirstPostorderForward)]
  #[case::breadth_first_forward(        Traversal::BreadthFirstForward)]
  #[case::breadth_first_backward(       Traversal::BreadthFirstBackward)]
  #[trace]
  fn test_graph_traverse_single_node(#[case] traversal: Traversal) -> Result<(), Report> {
    let tree = NamedGraph::new(&["root"], &[])?;

    let mut actual = vec![];
    traverse(&tree, traversal, |name| {
      actual.push(name);
      Ok(())
    })?;

    assert_eq!(vec!["root"], actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::depth_first_preorder_forward( Traversal::DepthFirstPreorderForward)]
  #[case::depth_first_postorder_forward(Traversal::DepthFirstPostorderForward)]
  #[case::breadth_first_forward(        Traversal::BreadthFirstForward)]
  #[case::breadth_first_backward(       Traversal::BreadthFirstBackward)]
  #[trace]
  fn test_graph_traverse_rejects_multiple_roots(#[case] traversal: Traversal) -> Result<(), Report> {
    let forest = NamedGraph::new(&["r1", "r2", "A", "B"], &[("r1", "A"), ("r2", "B")])?;

    let mut visited = vec![];
    let result = traverse(&forest, traversal, |name| {
      visited.push(name);
      Ok(())
    });

    assert_error!(
      result,
      "Graph must have exactly one root: Only trees with exactly one root are currently supported, but found 2 roots. This is an internal error. Please report it to developers."
    );
    assert_eq!(Vec::<&str>::new(), visited);
    Ok(())
  }

  mod helpers {
    use super::super::super::test_graph_support::tests::NamedGraph;
    use eyre::Report;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Traversal {
      DepthFirstPreorderForward,
      DepthFirstPostorderForward,
      BreadthFirstForward,
      BreadthFirstBackward,
    }

    pub(super) fn fixture_balanced() -> Result<NamedGraph, Report> {
      NamedGraph::new(
        &["root", "AB", "CD", "A", "B", "C", "D"],
        &[
          ("root", "AB"),
          ("root", "CD"),
          ("AB", "A"),
          ("AB", "B"),
          ("CD", "C"),
          ("CD", "D"),
        ],
      )
    }

    pub(super) fn traverse(
      tree: &NamedGraph,
      traversal: Traversal,
      mut visit: impl FnMut(&'static str) -> Result<(), Report>,
    ) -> Result<(), Report> {
      let graph = &tree.graph;
      match traversal {
        Traversal::DepthFirstPreorderForward => graph.iter_depth_first_preorder_forward(|node| visit(tree.name(node.key))),
        Traversal::DepthFirstPostorderForward => {
          graph.iter_depth_first_postorder_forward(|node| visit(tree.name(node.key)))
        },
        Traversal::BreadthFirstForward => graph.iter_breadth_first_forward(|node| visit(tree.name(node.key))),
        Traversal::BreadthFirstBackward => graph.iter_breadth_first_backward(|node| visit(tree.name(node.key))),
      }
    }
  }
}
