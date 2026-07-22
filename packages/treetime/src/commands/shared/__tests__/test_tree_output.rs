#[cfg(test)]
mod tests {
  use crate::commands::shared::tree_output::{
    TreeOutputAdapter, TreeOutputData, TreeOutputEdge, TreeOutputMutation, TreeOutputNode, TreeOutputTrait,
    format_number,
  };
  use crate::payload::ancestral::GraphAncestral;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
  use treetime_graph::topology_order::TopologyOrderSpec;
  use treetime_io::auspice::{auspice_from_graph, auspice_write_str};
  use treetime_io::auspice_types::AuspiceTree;
  use treetime_io::nwk::{EdgeToNwk, NodeToNwk, nwk_read_str};
  use treetime_io::phyloxml::{phyloxml_from_graph, phyloxml_write_str};
  use treetime_io::usher_mat::usher_from_graph;
  use treetime_primitives::AsciiChar;

  #[test]
  fn test_tree_output_all_schema_adapters_use_final_graph_order() -> Result<(), Report> {
    let mut graph: GraphAncestral = nwk_read_str("((A,B,C)large,D)root;")?;
    TopologyOrderSpec::descendant_count(false).apply(&mut graph)?;

    let auspice = auspice_from_graph::<TreeOutputAdapter, _, _, _>(&graph)?;
    let auspice_children = auspice
      .tree
      .children
      .iter()
      .map(|child| child.name.as_str())
      .collect_vec();

    let phyloxml = phyloxml_from_graph::<TreeOutputAdapter, _, _, _>(&graph)?;
    let phyloxml_children = phyloxml.phylogeny[0]
      .clade
      .as_ref()
      .expect("PhyloXML must contain a root clade")
      .clade
      .iter()
      .map(|child| child.name.as_deref().unwrap_or_default())
      .collect_vec();

    let usher = usher_from_graph::<TreeOutputAdapter, _, _, _>(&graph)?;
    let usher_graph: GraphAncestral = nwk_read_str(&usher.newick)?;
    let usher_children = root_child_names(&usher_graph)?;

    let expected = vec!["D", "large"];
    assert_eq!(expected, auspice_children);
    assert_eq!(expected, phyloxml_children);
    assert_eq!(expected, usher_children);

    Ok(())
  }

  #[test]
  fn test_tree_output_adapter_preserves_supported_graph_semantics() -> Result<(), Report> {
    let graph = helpers::semantic_graph()?;

    let json = auspice_write_str::<TreeOutputAdapter, _, _, _>(&graph)?;
    let auspice: AuspiceTree = serde_json::from_str(&json)?;
    let child_a = auspice.tree.children.iter().find(|child| child.name == "A").unwrap();
    let child_b = auspice.tree.children.iter().find(|child| child.name == "B").unwrap();
    assert_eq!(Some(0.5), child_a.node_attrs.div);
    assert_eq!(vec!["A1T".to_owned()], child_a.branch_attrs.mutations["nuc"]);
    let num_date = child_a.node_attrs.num_date.as_ref().unwrap();
    assert_ulps_eq!(2020.5, num_date.value, max_ulps = 0);
    assert_eq!(Some([2020.2, 2020.8]), num_date.confidence);
    assert_eq!("Yes", child_b.node_attrs.bad_branch.as_ref().unwrap().value);
    assert_eq!("No", auspice.tree.node_attrs.bad_branch.as_ref().unwrap().value);
    let coloring_keys = auspice
      .data
      .meta
      .colorings
      .iter()
      .map(|coloring| coloring.key.as_str())
      .collect_vec();
    assert!(coloring_keys.contains(&"num_date"));
    assert!(coloring_keys.contains(&"bad_branch"));
    assert!(coloring_keys.contains(&"country"));
    assert!(coloring_keys.contains(&"gt"));

    let xml = phyloxml_write_str::<TreeOutputAdapter, _, _, _>(&graph)?;
    assert!(xml.contains("<name>A</name>"));
    assert!(xml.contains("<branch_length>0.5</branch_length>"));
    assert!(xml.contains("treetime:divergence"));
    assert!(xml.contains("treetime:mutation"));

    let usher = usher_from_graph::<TreeOutputAdapter, _, _, _>(&graph)?;
    let mutations = usher
      .node_mutations
      .iter()
      .flat_map(|mutations| &mutations.mutation)
      .collect_vec();
    assert_eq!(1, mutations.len());
    assert_eq!(1, mutations[0].position);
    assert_eq!(0, mutations[0].ref_nuc);
    assert_eq!(0, mutations[0].par_nuc);
    assert_eq!(vec![3], mutations[0].mut_nuc);

    Ok(())
  }

  // Oracle: augur export_v2.format_number keeps `precision` significant figures in the
  // fractional part while preserving integer digits.
  #[test]
  fn test_tree_output_format_number_fractional_precision() {
    assert_ulps_eq!(0.123457, format_number(0.12345678, 6), max_ulps = 0);
    assert_ulps_eq!(123.456789, format_number(123.456789, 6), max_ulps = 0);
    assert_ulps_eq!(0.0, format_number(0.0, 6), max_ulps = 0);
    assert_ulps_eq!(2020.123, format_number(2020.1234567, 3), max_ulps = 0);
  }

  fn root_child_names<N, E, D>(graph: &Graph<N, E, D>) -> Result<Vec<String>, Report>
  where
    N: GraphNode + Named,
    E: GraphEdge,
    D: Send + Sync,
  {
    let root = graph.get_exactly_one_root()?;
    Ok(
      graph
        .children_of(&root.read_arc())
        .into_iter()
        .map(|(child, _)| {
          child
            .read_arc()
            .payload()
            .read_arc()
            .name()
            .expect("fixture child must have a name")
            .as_ref()
            .to_owned()
        })
        .collect(),
    )
  }

  mod helpers {
    use super::*;

    #[derive(Debug)]
    pub struct TestNode {
      name: String,
      div: f64,
      date: Option<f64>,
      bad_branch: bool,
    }

    impl GraphNode for TestNode {}

    impl Named for TestNode {
      fn name(&self) -> Option<impl AsRef<str>> {
        Some(&self.name)
      }

      fn set_name(&mut self, name: Option<impl AsRef<str>>) {
        self.name = name.map_or_else(String::new, |name| name.as_ref().to_owned());
      }
    }

    impl NodeToNwk for TestNode {
      fn nwk_name(&self) -> Option<impl AsRef<str>> {
        Some(&self.name)
      }
    }

    impl TreeOutputNode for TestNode {
      fn output_divergence(&self) -> Option<f64> {
        Some(self.div)
      }

      fn output_date(&self) -> Option<f64> {
        self.date
      }

      fn output_bad_branch(&self) -> bool {
        self.bad_branch
      }
    }

    #[derive(Debug)]
    pub struct TestEdge {
      branch_length: f64,
    }

    impl GraphEdge for TestEdge {}

    impl HasBranchLength for TestEdge {
      fn branch_length(&self) -> Option<f64> {
        Some(self.branch_length)
      }

      fn set_branch_length(&mut self, branch_length: Option<f64>) {
        if let Some(branch_length) = branch_length {
          self.branch_length = branch_length;
        }
      }
    }

    impl EdgeToNwk for TestEdge {
      fn nwk_weight(&self) -> Option<f64> {
        Some(self.branch_length)
      }
    }

    impl TreeOutputEdge for TestEdge {}

    #[derive(Debug, Default)]
    pub struct TestData {
      date_confidence: BTreeMap<GraphNodeKey, [f64; 2]>,
      mutations: BTreeMap<GraphEdgeKey, Vec<TreeOutputMutation>>,
      traits: BTreeMap<GraphNodeKey, BTreeMap<String, TreeOutputTrait>>,
    }

    impl TreeOutputData<TestNode, TestEdge> for TestData {
      fn trait_attributes(&self) -> Vec<String> {
        vec!["country".to_owned()]
      }

      fn has_dates(&self) -> bool {
        true
      }

      fn has_bad_branch(&self) -> bool {
        true
      }

      fn has_mutations(&self) -> bool {
        true
      }

      fn root_sequences(&self, _graph: &Graph<TestNode, TestEdge, Self>) -> Result<BTreeMap<String, String>, Report> {
        Ok(BTreeMap::from([("nuc".to_owned(), "A".to_owned())]))
      }

      fn date_confidence(&self, key: GraphNodeKey) -> Option<[f64; 2]> {
        self.date_confidence.get(&key).copied()
      }

      fn traits(&self, key: GraphNodeKey) -> BTreeMap<String, TreeOutputTrait> {
        self.traits.get(&key).cloned().unwrap_or_default()
      }

      fn mutations(
        &self,
        _graph: &Graph<TestNode, TestEdge, Self>,
        key: GraphEdgeKey,
      ) -> Result<Vec<TreeOutputMutation>, Report> {
        Ok(self.mutations.get(&key).cloned().unwrap_or_default())
      }
    }

    pub fn semantic_graph() -> Result<Graph<TestNode, TestEdge, TestData>, Report> {
      let mut graph = Graph::with_data(TestData::default());
      let root = graph.add_node(node("root", 0.0, Some(2020.0), false));
      let a = graph.add_node(node("A", 0.5, Some(2020.5), false));
      let b = graph.add_node(node("B", 1.0, None, true));
      let edge_a = graph.add_edge(root, a, TestEdge { branch_length: 0.5 })?;
      graph.add_edge(root, b, TestEdge { branch_length: 1.0 })?;
      graph.build()?;

      graph.data_mut().date_confidence.insert(a, [2020.2, 2020.8]);
      graph.data_mut().mutations.insert(
        edge_a,
        vec![TreeOutputMutation {
          gene: "nuc".to_owned(),
          position: 0,
          parent: AsciiChar::try_new(b'A')?,
          child: AsciiChar::try_new(b'T')?,
        }],
      );
      graph.data_mut().traits.insert(
        a,
        BTreeMap::from([(
          "country".to_owned(),
          TreeOutputTrait {
            value: "CH".to_owned(),
            ..TreeOutputTrait::default()
          },
        )]),
      );
      Ok(graph)
    }

    fn node(name: &str, div: f64, date: Option<f64>, bad_branch: bool) -> TestNode {
      TestNode {
        name: name.to_owned(),
        div,
        date,
        bad_branch,
      }
    }
  }
}
