#[cfg(test)]
mod tests {
  use crate::auspice::{auspice_read_str, auspice_write_str};
  use crate::auspice_types::AuspiceTree;
  use crate::tree_ir::auspice::{TreeIrAuspiceReader, TreeIrAuspiceWriter, format_number};
  use approx::assert_ulps_eq;
  use eyre::Report;
  use pretty_assertions::assert_eq;

  // Oracle: augur export_v2.format_number keeps `precision` significant figures in the
  // fractional part while preserving integer digits.
  #[test]
  fn test_tree_ir_auspice_format_number_fractional_precision() {
    assert_ulps_eq!(0.123457, format_number(0.12345678, 6), max_ulps = 0);
    assert_ulps_eq!(123.456789, format_number(123.456789, 6), max_ulps = 0);
    assert_ulps_eq!(0.0, format_number(0.0, 6), max_ulps = 0);
    assert_ulps_eq!(2020.123, format_number(2020.1234567, 3), max_ulps = 0);
  }

  #[test]
  fn test_tree_ir_auspice_div_is_cumulative_and_rounded() -> Result<(), Report> {
    let graph = helpers::sample_graph()?;
    let json = auspice_write_str::<TreeIrAuspiceWriter, _, _, _>(&graph)?;
    let tree: AuspiceTree = serde_json::from_str(&json)?;

    assert_eq!("root", tree.tree.name);
    assert_eq!(Some(0.0), tree.tree.node_attrs.div);
    let child_a = tree.tree.children.iter().find(|c| c.name == "A").unwrap();
    assert_eq!(Some(0.5), child_a.node_attrs.div);
    Ok(())
  }

  #[test]
  fn test_tree_ir_auspice_mutations_grouped_by_gene() -> Result<(), Report> {
    let graph = helpers::sample_graph()?;
    let json = auspice_write_str::<TreeIrAuspiceWriter, _, _, _>(&graph)?;
    let tree: AuspiceTree = serde_json::from_str(&json)?;

    let child_a = tree.tree.children.iter().find(|c| c.name == "A").unwrap();
    assert_eq!(vec!["A123T".to_owned()], child_a.branch_attrs.mutations["nuc"]);
    Ok(())
  }

  #[test]
  fn test_tree_ir_auspice_num_date_with_confidence() -> Result<(), Report> {
    let graph = helpers::sample_graph()?;
    let json = auspice_write_str::<TreeIrAuspiceWriter, _, _, _>(&graph)?;
    let tree: AuspiceTree = serde_json::from_str(&json)?;

    let child_a = tree.tree.children.iter().find(|c| c.name == "A").unwrap();
    let num_date = child_a.node_attrs.num_date.as_ref().unwrap();
    assert_ulps_eq!(2020.5, num_date.value, max_ulps = 0);
    assert_eq!(Some([2020.2, 2020.8]), num_date.confidence);
    Ok(())
  }

  #[test]
  fn test_tree_ir_auspice_bad_branch_categorical() -> Result<(), Report> {
    let graph = helpers::sample_graph()?;
    let json = auspice_write_str::<TreeIrAuspiceWriter, _, _, _>(&graph)?;
    let tree: AuspiceTree = serde_json::from_str(&json)?;

    let child_b = tree.tree.children.iter().find(|c| c.name == "B").unwrap();
    assert_eq!("Yes", child_b.node_attrs.bad_branch.as_ref().unwrap().value);
    assert_eq!("No", tree.tree.node_attrs.bad_branch.as_ref().unwrap().value);
    Ok(())
  }

  #[test]
  fn test_tree_ir_auspice_colorings_present() -> Result<(), Report> {
    let graph = helpers::sample_graph()?;
    let json = auspice_write_str::<TreeIrAuspiceWriter, _, _, _>(&graph)?;
    let tree: AuspiceTree = serde_json::from_str(&json)?;

    let keys: Vec<&str> = tree.data.meta.colorings.iter().map(|c| c.key.as_str()).collect();
    assert!(keys.contains(&"num_date"));
    assert!(keys.contains(&"bad_branch"));
    assert!(keys.contains(&"gt"));
    Ok(())
  }

  #[test]
  fn test_tree_ir_auspice_roundtrip_is_idempotent() -> Result<(), Report> {
    let graph = helpers::sample_graph()?;
    let json1 = auspice_write_str::<TreeIrAuspiceWriter, _, _, _>(&graph)?;
    let graph2 = auspice_read_str::<TreeIrAuspiceReader, _, _, _>(&json1)?;
    let json2 = auspice_write_str::<TreeIrAuspiceWriter, _, _, _>(&graph2)?;
    assert_eq!(json1, json2);
    Ok(())
  }

  #[test]
  fn test_tree_ir_auspice_non_finite_div_errors() -> Result<(), Report> {
    let graph = helpers::graph_with_root_div(f64::NAN)?;
    let err = auspice_write_str::<TreeIrAuspiceWriter, _, _, _>(&graph).unwrap_err();
    assert!(format!("{err:?}").contains("non-finite div"));
    Ok(())
  }

  mod helpers {
    use crate::graph::TreeIrGraph;
    use crate::tree_ir::mutation::{NUC_GENE, TreeIrSub};
    use crate::tree_ir::types::{TreeIrData, TreeIrEdge, TreeIrNode};
    use eyre::Report;
    use treetime_primitives::AsciiChar;

    pub fn nuc(c: char) -> AsciiChar {
      AsciiChar::try_from_char(c).unwrap()
    }

    pub fn sub(position: usize, parent: char, child: char) -> TreeIrSub {
      TreeIrSub {
        gene: NUC_GENE.to_owned(),
        position,
        parent: nuc(parent),
        child: nuc(child),
      }
    }

    pub fn data() -> TreeIrData {
      TreeIrData {
        has_dates: true,
        has_bad_branch: true,
        has_mutations: true,
        ..TreeIrData::default()
      }
    }

    pub fn sample_graph() -> Result<TreeIrGraph, Report> {
      let mut graph = TreeIrGraph::with_data(data());
      let root = graph.add_node(TreeIrNode {
        name: Some("root".to_owned()),
        div: Some(0.0),
        date: Some(2020.0),
        ..TreeIrNode::default()
      });
      let a = graph.add_node(TreeIrNode {
        name: Some("A".to_owned()),
        div: Some(0.5),
        date: Some(2020.5),
        date_confidence: Some([2020.2, 2020.8]),
        ..TreeIrNode::default()
      });
      let b = graph.add_node(TreeIrNode {
        name: Some("B".to_owned()),
        div: Some(1.0),
        bad_branch: true,
        ..TreeIrNode::default()
      });
      graph.add_edge(root, a, TreeIrEdge {
        branch_length: Some(0.5),
        mutations: vec![sub(123, 'A', 'T')],
        ..TreeIrEdge::default()
      })?;
      graph.add_edge(root, b, TreeIrEdge {
        branch_length: Some(1.0),
        ..TreeIrEdge::default()
      })?;
      graph.build()?;
      Ok(graph)
    }

    pub fn graph_with_root_div(div: f64) -> Result<TreeIrGraph, Report> {
      let mut graph = TreeIrGraph::with_data(data());
      graph.add_node(TreeIrNode {
        name: Some("root".to_owned()),
        div: Some(div),
        ..TreeIrNode::default()
      });
      graph.build()?;
      Ok(graph)
    }
  }
}
