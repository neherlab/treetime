#[cfg(test)]
mod tests {
  use crate::tree_ir::usher::{TreeIrUsherReader, TreeIrUsherWriter};
  use crate::usher_mat::{usher_from_graph, usher_to_graph};
  use eyre::Report;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_tree_ir_usher_encodes_nuc_substitution() -> Result<(), Report> {
    let graph = helpers::sample_graph(helpers::nuc_sub())?;
    let tree = usher_from_graph::<TreeIrUsherWriter, _, _, _>(&graph)?;

    let muts: Vec<_> = tree.node_mutations.iter().flat_map(|ml| &ml.mutation).collect();
    assert_eq!(1, muts.len());
    assert_eq!(123, muts[0].position);
    assert_eq!(0, muts[0].par_nuc); // A
    assert_eq!(vec![3], muts[0].mut_nuc); // T
    assert_eq!(0, muts[0].ref_nuc); // no reference -> falls back to parent (A)
    Ok(())
  }

  #[test]
  fn test_tree_ir_usher_roundtrip_preserves_mutations_and_topology() -> Result<(), Report> {
    let graph = helpers::sample_graph(helpers::nuc_sub())?;
    let tree = usher_from_graph::<TreeIrUsherWriter, _, _, _>(&graph)?;
    let graph2 = usher_to_graph::<TreeIrUsherReader, _, _, _>(&tree)?;
    let tree2 = usher_from_graph::<TreeIrUsherWriter, _, _, _>(&graph2)?;

    assert_eq!(tree.newick, tree2.newick);
    assert_eq!(tree.node_mutations, tree2.node_mutations);
    Ok(())
  }

  #[test]
  fn test_tree_ir_usher_drops_amino_acid_mutations() -> Result<(), Report> {
    let graph = helpers::sample_graph(helpers::aa_sub())?;
    let tree = usher_from_graph::<TreeIrUsherWriter, _, _, _>(&graph)?;
    let total: usize = tree.node_mutations.iter().map(|ml| ml.mutation.len()).sum();
    assert_eq!(0, total);
    Ok(())
  }

  #[test]
  fn test_tree_ir_usher_drops_indels() -> Result<(), Report> {
    let graph = helpers::sample_graph_with_indel()?;
    let tree = usher_from_graph::<TreeIrUsherWriter, _, _, _>(&graph)?;
    let total: usize = tree.node_mutations.iter().map(|ml| ml.mutation.len()).sum();
    assert_eq!(0, total);
    Ok(())
  }

  mod helpers {
    use crate::graph::TreeIrGraph;
    use crate::tree_ir::mutation::{IndelKind, NUC_GENE, TreeIrIndel, TreeIrSub};
    use crate::tree_ir::types::{TreeIrData, TreeIrEdge, TreeIrNode};
    use eyre::Report;
    use treetime_primitives::AsciiChar;

    fn nuc(c: char) -> AsciiChar {
      AsciiChar::try_from_char(c).unwrap()
    }

    pub fn nuc_sub() -> TreeIrSub {
      TreeIrSub {
        gene: NUC_GENE.to_owned(),
        position: 123,
        parent: nuc('A'),
        child: nuc('T'),
      }
    }

    pub fn aa_sub() -> TreeIrSub {
      TreeIrSub {
        gene: "M".to_owned(),
        position: 1,
        parent: nuc('M'),
        child: nuc('I'),
      }
    }

    fn named_node(name: &str) -> TreeIrNode {
      TreeIrNode {
        name: Some(name.to_owned()),
        ..TreeIrNode::default()
      }
    }

    pub fn sample_graph(sub: TreeIrSub) -> Result<TreeIrGraph, Report> {
      let mut graph = TreeIrGraph::with_data(TreeIrData::default());
      let root = graph.add_node(named_node("root"));
      let a = graph.add_node(named_node("A"));
      let b = graph.add_node(named_node("B"));
      graph.add_edge(root, a, TreeIrEdge {
        branch_length: Some(0.5),
        mutations: vec![sub],
        ..TreeIrEdge::default()
      })?;
      graph.add_edge(root, b, TreeIrEdge {
        branch_length: Some(1.0),
        ..TreeIrEdge::default()
      })?;
      graph.build()?;
      Ok(graph)
    }

    pub fn sample_graph_with_indel() -> Result<TreeIrGraph, Report> {
      let mut graph = TreeIrGraph::with_data(TreeIrData::default());
      let root = graph.add_node(named_node("root"));
      let a = graph.add_node(named_node("A"));
      graph.add_edge(root, a, TreeIrEdge {
        branch_length: Some(0.5),
        indels: vec![TreeIrIndel {
          gene: NUC_GENE.to_owned(),
          kind: IndelKind::Deletion,
          start: 10,
          seq: vec![nuc('A'), nuc('C')],
        }],
        ..TreeIrEdge::default()
      })?;
      graph.build()?;
      Ok(graph)
    }
  }
}
