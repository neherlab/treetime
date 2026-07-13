#[cfg(test)]
mod tests {
  use crate::graph::TreeIrGraph;
  use crate::phyloxml::{phyloxml_read_str, phyloxml_write_str};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use treetime_graph::node::Named;

  #[test]
  fn test_tree_ir_phyloxml_emits_name_branch_length_and_div_property() -> Result<(), Report> {
    let graph = helpers::sample_graph()?;
    let xml = phyloxml_write_str(&graph)?;
    assert!(xml.contains("<name>A</name>"));
    assert!(xml.contains("<branch_length>0.5</branch_length>"));
    assert!(xml.contains("treetime:divergence"));
    assert!(xml.contains("treetime:mutation"));
    Ok(())
  }

  #[test]
  fn test_tree_ir_phyloxml_roundtrip_is_idempotent() -> Result<(), Report> {
    let graph = helpers::sample_graph()?;
    let xml1 = phyloxml_write_str(&graph)?;
    let graph2: TreeIrGraph = phyloxml_read_str(&xml1)?;
    let xml2 = phyloxml_write_str(&graph2)?;
    assert_eq!(xml1, xml2);
    Ok(())
  }

  #[test]
  fn test_tree_ir_phyloxml_escapes_xml_special_chars_in_name() -> Result<(), Report> {
    let graph = helpers::single_node_named("a<b&c\"d")?;
    let xml = phyloxml_write_str(&graph)?;
    // quick-xml escapes special characters in element text.
    assert!(xml.contains("a&lt;b&amp;c"));
    let graph2: TreeIrGraph = phyloxml_read_str(&xml)?;
    let root = graph2.get_exactly_one_root()?;
    let name = root
      .read_arc()
      .payload()
      .read_arc()
      .name()
      .map(|n| n.as_ref().to_owned());
    assert_eq!(Some("a<b&c\"d".to_owned()), name);
    Ok(())
  }

  mod helpers {
    use crate::graph::TreeIrGraph;
    use crate::tree_ir::mutation::{NUC_GENE, TreeIrSub};
    use crate::tree_ir::types::{TreeIrData, TreeIrEdge, TreeIrNode};
    use eyre::Report;
    use treetime_primitives::AsciiChar;

    fn nuc(c: char) -> AsciiChar {
      AsciiChar::try_from_char(c).unwrap()
    }

    pub fn sample_graph() -> Result<TreeIrGraph, Report> {
      let mut graph = TreeIrGraph::with_data(TreeIrData::default());
      let root = graph.add_node(TreeIrNode {
        name: Some("root".to_owned()),
        div: Some(0.0),
        ..TreeIrNode::default()
      });
      let a = graph.add_node(TreeIrNode {
        name: Some("A".to_owned()),
        div: Some(0.5),
        date: Some(2020.5),
        ..TreeIrNode::default()
      });
      graph.add_edge(
        root,
        a,
        TreeIrEdge {
          branch_length: Some(0.5),
          mutations: vec![TreeIrSub {
            gene: NUC_GENE.to_owned(),
            position: 7,
            parent: nuc('A'),
            child: nuc('G'),
          }],
          ..TreeIrEdge::default()
        },
      )?;
      graph.build()?;
      Ok(graph)
    }

    pub fn single_node_named(name: &str) -> Result<TreeIrGraph, Report> {
      let mut graph = TreeIrGraph::with_data(TreeIrData::default());
      graph.add_node(TreeIrNode {
        name: Some(name.to_owned()),
        ..TreeIrNode::default()
      });
      graph.build()?;
      Ok(graph)
    }
  }
}
