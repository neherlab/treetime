#[cfg(test)]
mod tests {
  use crate::commands::timetree::output::ir::build_timetree_ir;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_io::auspice::auspice_write_str;
  use treetime_io::auspice_types::AuspiceTree;
  use treetime_io::tree_ir::auspice::TreeIrAuspiceWriter;

  // Without per-edge mutation counts, divergence is taken directly from the node's
  // stored cumulative substitutions-per-site value.
  #[test]
  fn test_build_timetree_ir_div_from_node_value() -> Result<(), Report> {
    let graph = helpers::two_node_graph()?;
    let ir = build_timetree_ir(&graph, None, None, None)?;
    let tree: AuspiceTree = serde_json::from_str(&auspice_write_str::<TreeIrAuspiceWriter, _, _, _>(&ir)?)?;

    assert_eq!(Some(0.0), tree.tree.node_attrs.div);
    let child = &tree.tree.children[0];
    assert_eq!(Some(0.5), child.node_attrs.div);
    Ok(())
  }

  // With per-edge mutation counts, divergence accumulates the counts from the root.
  #[test]
  fn test_build_timetree_ir_div_accumulates_mutation_counts() -> Result<(), Report> {
    let (graph, edge_key) = helpers::two_node_graph_with_edge_key()?;
    let counts = btreemap! { edge_key => 3_usize };
    let ir = build_timetree_ir(&graph, None, Some(&counts), None)?;
    let tree: AuspiceTree = serde_json::from_str(&auspice_write_str::<TreeIrAuspiceWriter, _, _, _>(&ir)?)?;

    assert_eq!(Some(0.0), tree.tree.node_attrs.div);
    assert_eq!(Some(3.0), tree.tree.children[0].node_attrs.div);
    Ok(())
  }

  #[test]
  fn test_build_timetree_ir_maps_date_and_bad_branch() -> Result<(), Report> {
    let graph = helpers::two_node_graph()?;
    let ir = build_timetree_ir(&graph, None, None, None)?;
    let tree: AuspiceTree = serde_json::from_str(&auspice_write_str::<TreeIrAuspiceWriter, _, _, _>(&ir)?)?;

    assert_ulps_eq!(
      2020.0,
      tree.tree.node_attrs.num_date.as_ref().unwrap().value,
      max_ulps = 0
    );
    assert_eq!("No", tree.tree.node_attrs.bad_branch.as_ref().unwrap().value);
    assert_eq!(
      "Yes",
      tree.tree.children[0].node_attrs.bad_branch.as_ref().unwrap().value
    );
    Ok(())
  }

  #[test]
  fn test_build_timetree_ir_emits_timetree_colorings() -> Result<(), Report> {
    let graph = helpers::two_node_graph()?;
    let ir = build_timetree_ir(&graph, None, None, None)?;
    let tree: AuspiceTree = serde_json::from_str(&auspice_write_str::<TreeIrAuspiceWriter, _, _, _>(&ir)?)?;

    let keys: Vec<&str> = tree.data.meta.colorings.iter().map(|c| c.key.as_str()).collect();
    assert!(keys.contains(&"num_date"));
    assert!(keys.contains(&"bad_branch"));
    assert_eq!(Some("bad_branch".to_owned()), tree.data.meta.display_defaults.color_by);
    Ok(())
  }

  mod helpers {
    use crate::partition::timetree::GraphTimetree;
    use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
    use eyre::Report;
    use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};
    use treetime_graph::node::Named;

    pub fn node(name: &str, div: f64, time: Option<f64>, bad_branch: bool) -> NodeTimetree {
      let mut node = NodeTimetree {
        time,
        bad_branch,
        div,
        ..NodeTimetree::default()
      };
      node.set_name(Some(name));
      node
    }

    pub fn edge(branch_length: f64) -> EdgeTimetree {
      let mut edge = EdgeTimetree::default();
      edge.set_branch_length(Some(branch_length));
      edge
    }

    pub fn two_node_graph() -> Result<GraphTimetree, Report> {
      Ok(two_node_graph_with_edge_key()?.0)
    }

    pub fn two_node_graph_with_edge_key() -> Result<(GraphTimetree, GraphEdgeKey), Report> {
      let mut graph = GraphTimetree::new();
      let root = graph.add_node(node("root", 0.0, Some(2020.0), false));
      let child = graph.add_node(node("A", 0.5, Some(2020.5), true));
      let edge_key = graph.add_edge(root, child, edge(0.5))?;
      graph.build()?;
      Ok((graph, edge_key))
    }
  }
}
