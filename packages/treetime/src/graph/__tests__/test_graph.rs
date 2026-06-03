#[cfg(test)]
mod tests {
  use std::collections::BTreeSet;
  use std::sync::Arc;

  use eyre::Report;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use treetime_utils::make_error;

  use treetime_graph::breadth_first::GraphTraversalContinuation;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::Named;
  use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};

  use crate::test_utils::{TestEdge, TestNode, find_edge_key, find_node_key_by_name};

  #[test]
  fn test_graph_traversal_serial_depth_first_preorder_forward() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let mut actual = vec![];
    graph.iter_depth_first_preorder_forward(|node| {
      actual.push(node.payload.name().unwrap().as_ref().to_owned());
      Ok(())
    })?;

    assert_eq!(vec!["root", "AB", "A", "B", "CD", "C", "D"], actual);

    Ok(())
  }

  #[test]
  fn test_graph_traversal_serial_depth_first_postorder_forward() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let mut actual = vec![];
    graph.iter_depth_first_postorder_forward(|node| {
      actual.push(node.payload.name().unwrap().as_ref().to_owned());
      Ok(())
    })?;

    assert_eq!(vec!["A", "B", "AB", "C", "D", "CD", "root"], actual);

    Ok(())
  }

  #[test]
  fn test_graph_traversal_serial_breadth_first_forward() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let mut actual = vec![];
    graph.iter_breadth_first_forward(|node| {
      actual.push(node.payload.name().unwrap().as_ref().to_owned());
      Ok(())
    })?;

    assert_eq!(vec!["root", "AB", "CD", "A", "B", "C", "D"], actual);

    Ok(())
  }

  #[test]
  fn test_graph_traversal_serial_breadth_first_reverse() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let mut actual = vec![];
    graph.iter_breadth_first_backward(|node| {
      actual.push(node.payload.name().unwrap().as_ref().to_owned());
      Ok(())
    })?;

    assert_eq!(vec!["D", "C", "B", "A", "CD", "AB", "root"], actual);

    Ok(())
  }

  #[test]
  fn test_graph_traversal_parallel_breadth_first_forward() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let actual = Arc::new(RwLock::new(vec![]));
    graph.par_iter_breadth_first_forward(|node| {
      actual
        .write_arc()
        .push(node.payload.name().unwrap().as_ref().to_owned());
      Ok(GraphTraversalContinuation::Continue)
    })?;

    assert_eq!(&vec!["root", "AB", "CD", "A", "B", "C", "D"], &*actual.read());

    Ok(())
  }

  #[test]
  fn test_graph_traversal_parallel_breadth_first_backward() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let actual = Arc::new(RwLock::new(vec![]));
    graph.par_iter_breadth_first_backward(|node| {
      actual
        .write_arc()
        .push(node.payload.name().unwrap().as_ref().to_owned());
      Ok(GraphTraversalContinuation::Continue)
    })?;

    assert_eq!(&vec!["D", "C", "B", "A", "CD", "AB", "root"], &*actual.read());

    Ok(())
  }

  #[test]
  fn test_graph_traversal_try_parallel_breadth_first_forward_ok() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let actual = Arc::new(RwLock::new(vec![]));
    graph.par_iter_breadth_first_forward(|node| {
      actual
        .write_arc()
        .push(node.payload.name().unwrap().as_ref().to_owned());
      Ok(GraphTraversalContinuation::Continue)
    })?;

    assert_eq!(&vec!["root", "AB", "CD", "A", "B", "C", "D"], &*actual.read());

    Ok(())
  }

  #[test]
  fn test_graph_traversal_try_parallel_breadth_first_forward_error() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    // Only one node errors, so the captured error is deterministic.
    let result = graph.par_iter_breadth_first_forward(|node| {
      let name = node.payload.name().unwrap().as_ref().to_owned();
      if name == "B" {
        return make_error!("boom {name}");
      }
      Ok(GraphTraversalContinuation::Continue)
    });

    assert_eq!("boom B", result.unwrap_err().to_string());

    Ok(())
  }

  #[test]
  fn test_graph_traversal_try_parallel_breadth_first_backward_ok() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let actual = Arc::new(RwLock::new(vec![]));
    graph.par_iter_breadth_first_backward(|node| {
      actual
        .write_arc()
        .push(node.payload.name().unwrap().as_ref().to_owned());
      Ok(GraphTraversalContinuation::Continue)
    })?;

    assert_eq!(&vec!["D", "C", "B", "A", "CD", "AB", "root"], &*actual.read());

    Ok(())
  }

  #[test]
  fn test_graph_traversal_try_parallel_breadth_first_backward_error() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let result = graph.par_iter_breadth_first_backward(|node| {
      let name = node.payload.name().unwrap().as_ref().to_owned();
      if name == "B" {
        return make_error!("boom {name}");
      }
      Ok(GraphTraversalContinuation::Continue)
    });

    assert_eq!("boom B", result.unwrap_err().to_string());

    Ok(())
  }

  #[test]
  fn test_graph_traversal_try_serial_breadth_first_forward_ok() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    // Serial traversal captures mutable outer state directly, with no Arc/Mutex.
    let mut visited = vec![];
    graph.iter_breadth_first_forward(|node| {
      visited.push(node.payload.name().unwrap().as_ref().to_owned());
      Ok(())
    })?;

    assert_eq!(vec!["root", "AB", "CD", "A", "B", "C", "D"], visited);

    Ok(())
  }

  #[test]
  fn test_graph_traversal_try_serial_breadth_first_forward_stops_at_first_error() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    // "AB" and "CD" both error, but "AB" is visited first in breadth-first order. The traversal
    // must surface the "AB" error and stop before visiting "CD" or any deeper node.
    let mut visited = vec![];
    let result = graph.iter_breadth_first_forward(|node| {
      let name = node.payload.name().unwrap().as_ref().to_owned();
      visited.push(name.clone());
      if name == "AB" || name == "CD" {
        return make_error!("boom {name}");
      }
      Ok(())
    });

    assert_eq!("boom AB", result.unwrap_err().to_string());
    assert_eq!(vec!["root", "AB"], visited);

    Ok(())
  }

  #[test]
  fn test_graph_traversal_try_serial_breadth_first_forward_single_node() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("root:0.01;")?;

    let mut visited = vec![];
    graph.iter_breadth_first_forward(|node| {
      visited.push(node.payload.name().unwrap().as_ref().to_owned());
      Ok(())
    })?;

    assert_eq!(vec!["root"], visited);

    Ok(())
  }

  #[test]
  fn test_graph_traversal_try_serial_breadth_first_forward_error_at_root() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1)root:0.01;")?;

    let mut visited = vec![];
    let result = graph.iter_breadth_first_forward(|node| {
      let name = node.payload.name().unwrap().as_ref().to_owned();
      visited.push(name.clone());
      if name == "root" {
        return make_error!("boom {name}");
      }
      Ok(())
    });

    assert_eq!("boom root", result.unwrap_err().to_string());
    assert_eq!(vec!["root"], visited);

    Ok(())
  }

  #[test]
  fn test_graph_traversal_try_serial_breadth_first_backward_ok() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let mut visited = vec![];
    graph.iter_breadth_first_backward(|node| {
      visited.push(node.payload.name().unwrap().as_ref().to_owned());
      Ok(())
    })?;

    assert_eq!(vec!["D", "C", "B", "A", "CD", "AB", "root"], visited);

    Ok(())
  }

  #[test]
  fn test_graph_traversal_try_serial_breadth_first_backward_stops_at_first_error() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let mut visited = vec![];
    let result = graph.iter_breadth_first_backward(|node| {
      let name = node.payload.name().unwrap().as_ref().to_owned();
      visited.push(name.clone());
      if name == "AB" {
        return make_error!("boom {name}");
      }
      Ok(())
    });

    assert_eq!("boom AB", result.unwrap_err().to_string());
    assert_eq!(vec!["D", "C", "B", "A", "CD", "AB"], visited);

    Ok(())
  }

  #[test]
  fn test_graph_traversal_try_depth_first_postorder_forward_ok() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let mut visited = vec![];
    graph.iter_depth_first_postorder_forward(|node| {
      visited.push(node.payload.name().unwrap().as_ref().to_owned());
      Ok(())
    })?;

    assert_eq!(vec!["A", "B", "AB", "C", "D", "CD", "root"], visited);

    Ok(())
  }

  #[test]
  fn test_graph_traversal_try_depth_first_postorder_forward_stops_at_first_error() -> Result<(), Report> {
    let graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let mut visited = vec![];
    let result = graph.iter_depth_first_postorder_forward(|node| {
      let name = node.payload.name().unwrap().as_ref().to_owned();
      visited.push(name.clone());
      if name == "AB" {
        return make_error!("boom {name}");
      }
      Ok(())
    });

    assert_eq!("boom AB", result.unwrap_err().to_string());
    assert_eq!(vec!["A", "B", "AB"], visited);

    Ok(())
  }

  #[test]
  fn test_graph_collapse_edge_simple_chain() -> Result<(), Report> {
    let mut graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.2)internal:0.1)root;")?;

    let root_to_internal = find_edge_key(&graph, "root", "internal").unwrap();
    graph.collapse_edge(root_to_internal)?;

    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.2)root;");

    Ok(())
  }

  #[test]
  fn test_graph_collapse_edge_binary_tree() -> Result<(), Report> {
    let mut graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.2)internal:0.1,B:0.3)root;")?;

    let root_to_internal = find_edge_key(&graph, "root", "internal").unwrap();
    graph.collapse_edge(root_to_internal)?;

    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(B:0.3,A:0.2)root;");

    Ok(())
  }

  #[test]
  fn test_graph_collapse_edge_complex_tree() -> Result<(), Report> {
    let mut graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.3,B:0.4)left:0.1,(C:0.5,D:0.6)right:0.2)root;")?;

    let root_to_left = find_edge_key(&graph, "root", "left").unwrap();
    graph.collapse_edge(root_to_left)?;

    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "((C:0.5,D:0.6)right:0.2,A:0.3,B:0.4)root;");

    Ok(())
  }

  #[allow(clippy::assertions_on_result_states)]
  #[test]
  fn test_graph_collapse_edge_invalid_edge() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();

    let root_node = graph.add_node(TestNode(Some("root".to_owned())));
    let leaf_node = graph.add_node(TestNode(Some("A".to_owned())));
    graph.add_edge(root_node, leaf_node, TestEdge(Some(0.1)))?;
    graph.build()?;

    let invalid_key = GraphEdgeKey(9999);
    let result = graph.collapse_edge(invalid_key);

    assert!(result.is_err());

    Ok(())
  }

  #[test]
  fn test_graph_collapse_edge_leaf_edge() -> Result<(), Report> {
    let mut graph = nwk_read_str::<TestNode, TestEdge, ()>("(A:0.1,B:0.2)root;")?;

    let root_to_a = find_edge_key(&graph, "root", "A").unwrap();
    graph.collapse_edge(root_to_a)?;

    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(B:0.2)root;");

    Ok(())
  }

  #[test]
  fn test_graph_collapse_edge_no_duplicate_edges() -> Result<(), Report> {
    let mut graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.2)internal:0.1)root;")?;

    let root_to_internal = find_edge_key(&graph, "root", "internal").unwrap();
    graph.collapse_edge(root_to_internal)?;

    let root_key = find_node_key_by_name(&graph, "root").unwrap();
    let root_node = graph.get_node(root_key).unwrap();
    let root_node = root_node.read_arc();

    let outbound_edges = root_node.outbound();
    let unique_outbound: BTreeSet<_> = outbound_edges.iter().collect();
    assert_eq!(
      outbound_edges.len(),
      unique_outbound.len(),
      "Duplicate outbound edges detected"
    );

    let inbound_edges = root_node.inbound();
    let unique_inbound: BTreeSet<_> = inbound_edges.iter().collect();
    assert_eq!(
      inbound_edges.len(),
      unique_inbound.len(),
      "Duplicate inbound edges detected"
    );

    assert_eq!(outbound_edges.len(), 1);

    Ok(())
  }

  #[test]
  fn test_graph_collapse_edge_adjacency_lists_maintained() -> Result<(), Report> {
    let mut graph = nwk_read_str::<TestNode, TestEdge, ()>("((A:0.3,B:0.4)left:0.1,right:0.2)root;")?;

    let root_to_left = find_edge_key(&graph, "root", "left").unwrap();
    graph.collapse_edge(root_to_left)?;

    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(right:0.2,A:0.3,B:0.4)root;");

    Ok(())
  }

  #[test]
  fn test_graph_collapse_edge_multiple_inbound_edges() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdge, ()>::new();

    let source1_node = graph.add_node(TestNode(Some("source1".to_owned())));
    let source2_node = graph.add_node(TestNode(Some("source2".to_owned())));
    let target_node = graph.add_node(TestNode(Some("target".to_owned())));
    let leaf_node = graph.add_node(TestNode(Some("leaf".to_owned())));

    let edge_to_collapse = graph.add_edge(source1_node, target_node, TestEdge(Some(0.1)))?;
    let _other_inbound = graph.add_edge(source2_node, target_node, TestEdge(Some(0.2)))?;
    let _outbound_edge = graph.add_edge(target_node, leaf_node, TestEdge(Some(0.3)))?;

    graph.build()?;

    graph.collapse_edge(edge_to_collapse)?;

    let source1_ref = graph.get_node(source1_node).unwrap();
    let source1_binding = source1_ref.read_arc();
    let source1_outbound = source1_binding.outbound();
    assert_eq!(source1_outbound.len(), 1);

    if let Some(edge) = graph.get_edge(source1_outbound[0]) {
      let edge_ref = edge.read_arc();
      assert_eq!(edge_ref.source(), source1_node);
      assert_eq!(edge_ref.target(), leaf_node);
    }

    let source1_inbound = source1_binding.inbound();
    assert_eq!(source1_inbound.len(), 1);

    if let Some(edge) = graph.get_edge(source1_inbound[0]) {
      let edge_ref = edge.read_arc();
      assert_eq!(edge_ref.source(), source2_node);
      assert_eq!(edge_ref.target(), source1_node);
    }

    assert!(graph.get_node(target_node).is_none());

    Ok(())
  }

  #[test]
  fn test_graph_collapse_edge_adjacency_consistency() -> Result<(), Report> {
    let mut graph = nwk_read_str::<TestNode, TestEdge, ()>("((leaf1:0.2,leaf2:0.3)internal:0.1)root;")?;

    let root_to_internal = find_edge_key(&graph, "root", "internal").unwrap();
    graph.collapse_edge(root_to_internal)?;

    for node in graph.get_nodes() {
      let node_ref = node.read_arc();

      for &edge_key in node_ref.outbound() {
        let edge = graph.get_edge(edge_key).expect("Outbound edge should exist");
        let edge_ref = edge.read_arc();
        assert_eq!(edge_ref.source(), node_ref.key(), "Edge source should match node");
      }

      for &edge_key in node_ref.inbound() {
        let edge = graph.get_edge(edge_key).expect("Inbound edge should exist");
        let edge_ref = edge.read_arc();
        assert_eq!(edge_ref.target(), node_ref.key(), "Edge target should match node");
      }
    }

    Ok(())
  }
}
