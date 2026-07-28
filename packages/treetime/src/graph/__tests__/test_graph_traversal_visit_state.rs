#[cfg(test)]
mod tests {
  use std::collections::BTreeSet;
  use std::sync::{Arc, Barrier};
  use std::thread;

  use eyre::Report;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use rayon::ThreadPoolBuilder;
  use treetime_graph::breadth_first::{GraphTraversalContinuation, directed_breadth_first_traversal_forward};
  use treetime_graph::find_paths::exists_forward_path_between;
  use treetime_graph::node::GraphNodeKey;
  use treetime_utils::{assert_error, make_error};

  use self::helpers::graph_with_edges;

  #[test]
  fn test_graph_path_queries_are_sequentially_independent() -> Result<(), Report> {
    let (graph, nodes) = graph_with_edges(3, &[(0, 1), (1, 2)])?;

    let root_to_leaf = exists_forward_path_between(&graph, &nodes[0], &nodes[2]);
    let root_to_internal = exists_forward_path_between(&graph, &nodes[0], &nodes[1]);
    let internal_to_leaf = exists_forward_path_between(&graph, &nodes[1], &nodes[2]);

    // Reachability follows the directed chain regardless of prior queries.
    assert_eq!((true, true, true), (root_to_leaf, root_to_internal, internal_to_leaf));

    Ok(())
  }

  #[test]
  fn test_graph_traversals_are_concurrently_independent() -> Result<(), Report> {
    let (graph, _nodes) = graph_with_edges(3, &[(0, 1), (1, 2)])?;
    let pool = ThreadPoolBuilder::new().num_threads(2).build()?;
    let barrier = Arc::new(Barrier::new(2));
    let forward = Arc::new(RwLock::new(vec![]));
    let backward = Arc::new(RwLock::new(vec![]));

    let (forward_result, backward_result) = thread::scope(|scope| {
      let forward_result = scope.spawn(|| {
        pool.install(|| {
          graph.par_iter_breadth_first_forward(|node| {
            if node.is_root {
              barrier.wait();
            }
            forward.write().push(node.key);
            Ok(GraphTraversalContinuation::Continue)
          })
        })
      });
      let backward_result = scope.spawn(|| {
        pool.install(|| {
          graph.par_iter_breadth_first_backward(|node| {
            if node.is_leaf {
              barrier.wait();
            }
            backward.write().push(node.key);
            Ok(GraphTraversalContinuation::Continue)
          })
        })
      });
      (forward_result.join().unwrap(), backward_result.join().unwrap())
    });
    forward_result?;
    backward_result?;

    assert_eq!(vec![GraphNodeKey(0), GraphNodeKey(1), GraphNodeKey(2)], *forward.read());
    assert_eq!(
      vec![GraphNodeKey(2), GraphNodeKey(1), GraphNodeKey(0)],
      *backward.read()
    );

    Ok(())
  }

  #[test]
  fn test_graph_traversal_error_does_not_affect_later_query() -> Result<(), Report> {
    let (graph, nodes) = graph_with_edges(3, &[(0, 1), (1, 2)])?;
    let middle = nodes[1].read_arc().key();

    let result = directed_breadth_first_traversal_forward(&graph, &[Arc::clone(&nodes[0])], |node| {
      if node.key() == middle {
        return make_error!("injected traversal failure");
      }
      Ok(GraphTraversalContinuation::Continue)
    });
    assert_error!(result, "injected traversal failure");

    let actual = exists_forward_path_between(&graph, &nodes[0], &nodes[2]);
    assert!(actual);

    Ok(())
  }

  #[test]
  fn test_graph_traversal_stop_does_not_affect_later_traversal() -> Result<(), Report> {
    let (graph, nodes) = graph_with_edges(3, &[(0, 1), (1, 2)])?;
    let middle = nodes[1].read_arc().key();
    directed_breadth_first_traversal_forward(&graph, &[Arc::clone(&nodes[0])], |node| {
      Ok(if node.key() == middle {
        GraphTraversalContinuation::Stop
      } else {
        GraphTraversalContinuation::Continue
      })
    })?;

    let actual = Arc::new(RwLock::new(vec![]));
    directed_breadth_first_traversal_forward(&graph, &[Arc::clone(&nodes[0])], |node| {
      actual.write().push(node.key());
      Ok(GraphTraversalContinuation::Continue)
    })?;

    assert_eq!(vec![GraphNodeKey(0), GraphNodeKey(1), GraphNodeKey(2)], *actual.read());

    Ok(())
  }

  #[test]
  fn test_graph_traversal_visits_every_forest_node() -> Result<(), Report> {
    let (graph, _nodes) = graph_with_edges(4, &[(0, 1), (2, 3)])?;
    let actual = Arc::new(RwLock::new(BTreeSet::new()));

    graph.par_iter_breadth_first_forward(|node| {
      actual.write().insert(node.key);
      Ok(GraphTraversalContinuation::Continue)
    })?;

    let expected = BTreeSet::from([GraphNodeKey(0), GraphNodeKey(1), GraphNodeKey(2), GraphNodeKey(3)]);
    assert_eq!(expected, *actual.read());

    Ok(())
  }

  #[test]
  fn test_graph_traversal_terminates_on_cycle() -> Result<(), Report> {
    let (graph, nodes) = graph_with_edges(3, &[(0, 1), (1, 2), (2, 0)])?;
    let actual = Arc::new(RwLock::new(vec![]));

    directed_breadth_first_traversal_forward(&graph, &[Arc::clone(&nodes[0])], |node| {
      actual.write().push(node.key());
      Ok(GraphTraversalContinuation::Continue)
    })?;

    assert_eq!(vec![GraphNodeKey(0), GraphNodeKey(1), GraphNodeKey(2)], *actual.read());

    Ok(())
  }

  mod helpers {
    use eyre::Report;
    use treetime_graph::graph::{Graph, SafeNode};

    use crate::test_utils::{TestEdge, TestNode};

    pub fn graph_with_edges(
      node_count: usize,
      edges: &[(usize, usize)],
    ) -> Result<(Graph<TestNode, TestEdge, ()>, Vec<SafeNode<TestNode>>), Report> {
      let mut graph = Graph::new();
      let keys = (0..node_count)
        .map(|index| graph.add_node(TestNode(Some(format!("node_{index}")))))
        .collect::<Vec<_>>();
      edges
        .iter()
        .try_for_each(|(source, target)| graph.add_edge(keys[*source], keys[*target], TestEdge(None)).map(|_| ()))?;
      graph.build()?;
      let nodes = keys.iter().map(|key| graph.get_node(*key).unwrap()).collect::<Vec<_>>();
      Ok((graph, nodes))
    }
  }
}
