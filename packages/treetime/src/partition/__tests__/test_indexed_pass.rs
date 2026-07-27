#[cfg(test)]
mod tests {
  use crate::partition::indexed_pass::{IndexedPass, with_indexed_graph_payloads};
  use crate::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::{assert_error, make_report, o};

  use self::helpers::{pass_values, values_by_name};

  #[test]
  fn test_indexed_pass_backward_visits_children_before_parent() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("((A:1,B:1)AB:1,C:1)root;")?;
    let (mut nodes, mut edges) = pass_values(&graph);
    let mut pass = IndexedPass::new(&graph, &mut nodes, &mut edges, |_| Ok(0))?;

    pass.try_for_each_backward(|dependencies, slot| {
      let graph_node = graph.get_node(slot.key).expect("Indexed node must exist");
      slot.node = graph
        .children_of(&graph_node.read_arc())
        .iter()
        .map(|(child, _)| dependencies.node(child.read_arc().key()))
        .sum::<usize>()
        + 1;
      Ok(())
    })?;
    let (nodes, _) = pass.into_maps()?;
    let actual = values_by_name(&graph, &nodes);
    let expected = btreemap! {
      o!("A") => 1,
      o!("AB") => 3,
      o!("B") => 1,
      o!("C") => 1,
      o!("root") => 5,
    };

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_indexed_pass_forward_visits_parent_before_children() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("((A:1,B:1)AB:1,C:1)root;")?;
    let (mut nodes, mut edges) = pass_values(&graph);
    let mut pass = IndexedPass::new(&graph, &mut nodes, &mut edges, |_| Ok(0))?;

    pass.try_for_each_forward(|dependencies, slot| {
      slot.node = slot.parent_key.map_or(0, |parent| dependencies.node(parent) + 1);
      Ok(())
    })?;
    let (nodes, _) = pass.into_maps()?;
    let actual = values_by_name(&graph, &nodes);
    let expected = btreemap! {
      o!("A") => 2,
      o!("AB") => 1,
      o!("B") => 2,
      o!("C") => 1,
      o!("root") => 0,
    };

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_indexed_pass_roundtrip_preserves_all_values() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("((A:1,B:2)AB:3,C:4)root;")?;
    let mut nodes = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let key = node.read_arc().key();
        (key, key.as_usize())
      })
      .collect::<BTreeMap<_, _>>();
    let mut edges = graph
      .get_edges()
      .iter()
      .map(|edge| {
        let key = edge.read_arc().key();
        (key, key.as_usize())
      })
      .collect::<BTreeMap<_, _>>();
    let expected_nodes = nodes.clone();
    let expected_edges = edges.clone();

    let pass = IndexedPass::new(&graph, &mut nodes, &mut edges, |_| {
      unreachable!("all graph nodes are present")
    })?;
    let (actual_nodes, actual_edges) = pass.into_maps()?;

    assert_eq!(expected_nodes, actual_nodes);
    assert_eq!(expected_edges, actual_edges);
    Ok(())
  }

  #[test]
  fn test_indexed_pass_error_restores_all_values() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("((A:1,B:2)AB:3,C:4)root;")?;
    let mut nodes = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let key = node.read_arc().key();
        (key, key.as_usize())
      })
      .collect::<BTreeMap<_, _>>();
    let mut edges = graph
      .get_edges()
      .iter()
      .map(|edge| {
        let key = edge.read_arc().key();
        (key, key.as_usize())
      })
      .collect::<BTreeMap<_, _>>();
    let expected_nodes = nodes.clone();
    let expected_edges = edges.clone();
    let mut pass = IndexedPass::new(&graph, &mut nodes, &mut edges, |_| {
      unreachable!("all graph nodes are present")
    })?;

    let result = pass.try_for_each_backward(|_, _| Err(make_report!("injected pass failure")));
    assert_error!(result, "injected pass failure");
    let (actual_nodes, actual_edges) = pass.into_maps()?;

    assert_eq!(expected_nodes, actual_nodes);
    assert_eq!(expected_edges, actual_edges);
    Ok(())
  }

  #[test]
  fn test_indexed_graph_payloads_error_restores_graph() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("((A:1,B:2)AB:3,C:4)root;")?;
    let expected_nodes = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let node = node.read_arc();
        (node.key(), node.payload().read_arc().name.clone())
      })
      .collect::<BTreeMap<_, _>>();
    let expected_edges = graph
      .get_edges()
      .iter()
      .map(|edge| {
        let edge = edge.read_arc();
        (edge.key(), edge.payload().read_arc().branch_length)
      })
      .collect::<BTreeMap<_, _>>();

    let result = with_indexed_graph_payloads(&graph, |_| Err::<(), _>(make_report!("injected graph pass failure")));
    assert_error!(result, "injected graph pass failure");
    let actual_nodes = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let node = node.read_arc();
        (node.key(), node.payload().read_arc().name.clone())
      })
      .collect::<BTreeMap<_, _>>();
    let actual_edges = graph
      .get_edges()
      .iter()
      .map(|edge| {
        let edge = edge.read_arc();
        (edge.key(), edge.payload().read_arc().branch_length)
      })
      .collect::<BTreeMap<_, _>>();

    assert_eq!(expected_nodes, actual_nodes);
    assert_eq!(expected_edges, actual_edges);
    Ok(())
  }

  mod helpers {
    use crate::payload::ancestral::GraphAncestral;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::node::GraphNodeKey;

    pub fn pass_values(graph: &GraphAncestral) -> (BTreeMap<GraphNodeKey, usize>, BTreeMap<GraphEdgeKey, usize>) {
      let nodes = graph
        .get_nodes()
        .iter()
        .map(|node| (node.read_arc().key(), 0))
        .collect();
      let edges = graph
        .get_edges()
        .iter()
        .map(|edge| (edge.read_arc().key(), 0))
        .collect();
      (nodes, edges)
    }

    pub fn values_by_name(graph: &GraphAncestral, values: &BTreeMap<GraphNodeKey, usize>) -> BTreeMap<String, usize> {
      graph
        .get_nodes()
        .iter()
        .map(|node| {
          let node = node.read_arc();
          let name = node
            .payload()
            .read_arc()
            .name
            .as_ref()
            .expect("Fixture nodes must be named")
            .clone();
          (name, values[&node.key()])
        })
        .collect()
    }
  }
}
