#[cfg(test)]
mod tests {
  use crate::partition::indexed_pass::IndexedPass;
  use crate::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::{assert_error, make_report};

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

    let result = pass.try_for_each_backward_frontier(|_, _, _, _, _| Err(make_report!("injected pass failure")));
    assert_error!(result, "injected pass failure");
    let (actual_nodes, actual_edges) = pass.into_maps()?;

    assert_eq!(expected_nodes, actual_nodes);
    assert_eq!(expected_edges, actual_edges);
    Ok(())
  }
}
