use crate::convert::convert::{ConverterData, ConverterEdge, ConverterNode};
use crate::convert::mutation::{PartitionedMutations, format_mutation_list, parse_mutation_list};
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use serde_json::Value;
use std::collections::BTreeMap;
use treetime::graph::graph::Graph;
use treetime::graph::node::GraphNodeKey;
use treetime::io::auspice::{
  AuspiceGraphContext, AuspiceRead, AuspiceTree, AuspiceTreeBranchAttrs, AuspiceTreeContext, AuspiceTreeData,
  AuspiceTreeNode, AuspiceTreeNodeAttrs, AuspiceWrite,
};
use treetime::make_internal_report;

fn parse_auspice_mutations(auspice_mutations: &BTreeMap<String, Vec<String>>) -> Result<PartitionedMutations, Report> {
  auspice_mutations
    .iter()
    .map(|(partition, muts)| Ok((partition.clone(), parse_mutation_list(muts)?)))
    .try_collect()
}

fn format_auspice_mutations(mutations: &PartitionedMutations) -> BTreeMap<String, Vec<String>> {
  mutations
    .iter()
    .map(|(partition, muts)| (partition.clone(), format_mutation_list(muts)))
    .collect()
}

pub struct AuspiceWriter {
  pub divs: BTreeMap<GraphNodeKey, f64>,
}

impl AuspiceWrite<ConverterNode, ConverterEdge, ConverterData> for AuspiceWriter {
  fn new(_: &Graph<ConverterNode, ConverterEdge, ConverterData>) -> Result<Self, Report> {
    Ok(Self { divs: btreemap! {} })
  }

  fn auspice_data_from_graph_data(
    &self,
    graph: &Graph<ConverterNode, ConverterEdge, ConverterData>,
  ) -> Result<AuspiceTreeData, Report> {
    let data = graph.data().read_arc();
    Ok(AuspiceTreeData {
      version: data.version.clone(),
      meta: data.meta.clone(),
      root_sequence: data.root_sequence.clone(),
      other: data.other.clone(),
    })
  }

  fn auspice_node_from_graph_components(
    &mut self,
    context: &AuspiceGraphContext<ConverterNode, ConverterEdge, ConverterData>,
  ) -> Result<AuspiceTreeNode, Report> {
    let &AuspiceGraphContext {
      node_key,
      node,
      parent_key,
      parent,
      edge,
      graph,
    } = context;

    let name = node
      .name
      .as_ref()
      .ok_or_else(|| make_internal_report!("Encountered node with empty name"))?
      .to_owned();

    let parent_div = parent_key
      .and_then(|parent_key| self.divs.get(&parent_key).copied())
      .or(Some(0.0)); // FIXME: we probably should not assume 0 if parent div is missing or if it's a root node

    let branch_length = edge.and_then(|edge| edge.weight);

    let div = if let (Some(parent_div), Some(branch_length)) = (parent_div, branch_length) {
      Some(parent_div + branch_length)
    } else {
      Some(0.0) // FIXME: we probably should not assume 0 if div or branch_length is missing
    };

    if let Some(div) = div {
      self.divs.insert(node_key, div);
    }

    let mutations = edge.map_or(BTreeMap::new(), |edge| format_auspice_mutations(&edge.mutations));

    Ok(AuspiceTreeNode {
      name,
      branch_attrs: AuspiceTreeBranchAttrs {
        mutations,
        labels: None,
        other: Value::default(),
      },
      node_attrs: AuspiceTreeNodeAttrs {
        div,
        clade_membership: None,
        region: None,
        country: None,
        division: None,
        other: Value::default(),
      },
      children: vec![],
      other: Value::default(),
    })
  }
}

pub struct AuspiceReader {
  has_mutations: bool,
}

impl AuspiceRead<ConverterNode, ConverterEdge, ConverterData> for AuspiceReader {
  fn new(tree: &AuspiceTree) -> Result<Self, Report> {
    let has_mutations = tree
      .iter_depth_first_preorder()
      .any(|(_, node)| !node.branch_attrs.mutations.is_empty());
    Ok(Self { has_mutations })
  }

  fn auspice_data_to_graph_data(&mut self, tree: &AuspiceTree) -> Result<ConverterData, Report> {
    Ok(ConverterData {
      rooted: true,
      has_mutations: self.has_mutations,
      version: tree.data.version.clone(),
      meta: tree.data.meta.clone(),
      root_sequence: tree.data.root_sequence.clone(),
      other: tree.data.other.clone(),
    })
  }

  fn auspice_node_to_graph_components(
    &mut self,
    context: &AuspiceTreeContext,
  ) -> Result<(ConverterNode, ConverterEdge), Report> {
    let AuspiceTreeContext { node, .. } = context;
    let branch_length = context.branch_length();
    let mutations = parse_auspice_mutations(&node.branch_attrs.mutations)?;
    Ok((
      ConverterNode {
        name: Some(node.name.clone()),
      },
      ConverterEdge {
        weight: branch_length,
        mutations,
      },
    ))
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use treetime::io::auspice::{auspice_read_str, auspice_write_str};
  use treetime::io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};

  #[test]
  fn test_auspice_from_nwk() -> Result<(), Report> {
    let graph = nwk_read_str::<ConverterNode, ConverterEdge, ConverterData>(NWK_TEST_EXAMPLE)?;
    let output = auspice_write_str::<AuspiceWriter, _, _, _>(&graph)?;
    assert_eq!(AUSPICE_TEST_EXAMPLE, output);
    Ok(())
  }

  #[test]
  fn test_auspice_to_nwk() -> Result<(), Report> {
    let graph = auspice_read_str::<AuspiceReader, _, _, _>(AUSPICE_TEST_EXAMPLE)?;
    let output = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(NWK_TEST_EXAMPLE, output);
    Ok(())
  }

  #[test]
  fn test_auspice_roundtrip() -> Result<(), Report> {
    let input = indoc!(
      // language=json
      r#"{
        "version": "v2",
        "meta": {
          "description": "This is a test!",
          "name": "Test"
        },
        "root_sequence": {
          "nuc": "ACGTACGTACGTACGTACGTACGT"
        },
        "unknown_field": 42,
        "tree": {
          "name": "root",
          "node_attrs": {
            "div": 0.0
          },
          "children": [
            {
              "name": "AB",
              "node_attrs": {
                "div": 3.0
              },
              "children": [
                {
                  "name": "A",
                  "node_attrs": {
                    "div": 8.0
                  }
                },
                {
                  "name": "B",
                  "node_attrs": {
                    "div": 5.0
                  }
                }
              ]
            },
            {
              "name": "C",
              "node_attrs": {
                "div": 7.0
              }
            }
          ]
        }
      }"#
    );

    let graph = auspice_read_str::<AuspiceReader, _, _, _>(input)?;
    let output = auspice_write_str::<AuspiceWriter, _, _, _>(&graph)?;
    assert_eq!(input, output);
    Ok(())
  }

  const NWK_TEST_EXAMPLE: &str =
    "((A:0.8,(B:0.6,C:0.9,D:0.4)BCD:0.7)ABCD:0.5,(E:0.7,F:0.8)EF:0.9,(G:0.5,H:1.3,I:0.2,J:0.9)GHIJ:0.6)ABCDEFGHIJ;";

  const AUSPICE_TEST_EXAMPLE: &str = indoc!(
    // language=json
    r#"{
          "meta": {},
          "tree": {
            "name": "ABCDEFGHIJ",
            "node_attrs": {
              "div": 0.0
            },
            "children": [
              {
                "name": "ABCD",
                "node_attrs": {
                  "div": 0.5
                },
                "children": [
                  {
                    "name": "A",
                    "node_attrs": {
                      "div": 1.300000011920929
                    }
                  },
                  {
                    "name": "BCD",
                    "node_attrs": {
                      "div": 1.199999988079071
                    },
                    "children": [
                      {
                        "name": "B",
                        "node_attrs": {
                          "div": 1.800000011920929
                        }
                      },
                      {
                        "name": "C",
                        "node_attrs": {
                          "div": 2.099999964237213
                        }
                      },
                      {
                        "name": "D",
                        "node_attrs": {
                          "div": 1.5999999940395355
                        }
                      }
                    ]
                  }
                ]
              },
              {
                "name": "EF",
                "node_attrs": {
                  "div": 0.8999999761581421
                },
                "children": [
                  {
                    "name": "E",
                    "node_attrs": {
                      "div": 1.5999999642372131
                    }
                  },
                  {
                    "name": "F",
                    "node_attrs": {
                      "div": 1.699999988079071
                    }
                  }
                ]
              },
              {
                "name": "GHIJ",
                "node_attrs": {
                  "div": 0.6000000238418579
                },
                "children": [
                  {
                    "name": "G",
                    "node_attrs": {
                      "div": 1.100000023841858
                    }
                  },
                  {
                    "name": "H",
                    "node_attrs": {
                      "div": 1.899999976158142
                    }
                  },
                  {
                    "name": "I",
                    "node_attrs": {
                      "div": 0.8000000268220901
                    }
                  },
                  {
                    "name": "J",
                    "node_attrs": {
                      "div": 1.5
                    }
                  }
                ]
              }
            ]
          }
        }"#
  );

  #[test]
  fn test_auspice_mutation_roundtrip() -> Result<(), Report> {
    let input = indoc!(
      // language=json
      r#"{
        "meta": {},
        "tree": {
          "name": "root",
          "node_attrs": {
            "div": 0.0
          },
          "children": [
            {
              "name": "A",
              "branch_attrs": {
                "mutations": {
                  "nuc": ["A100T", "C200G"],
                  "S": ["D614G"]
                }
              },
              "node_attrs": {
                "div": 2.0
              }
            },
            {
              "name": "B",
              "branch_attrs": {
                "mutations": {
                  "nuc": ["G300A"]
                }
              },
              "node_attrs": {
                "div": 1.0
              }
            }
          ]
        }
      }"#
    );

    let graph = auspice_read_str::<AuspiceReader, _, _, _>(input)?;

    assert!(graph.data().read_arc().has_mutations);

    let edges = graph.get_edges();
    let edge_a = edges.iter().find(|e| {
      let edge = e.read_arc();
      let target_node = graph.get_node(edge.target());
      target_node.is_some_and(|n| n.read_arc().payload().read_arc().name.as_deref() == Some("A"))
    });
    let edge_a = edge_a.map(|e| e.read_arc());
    let edge_a_payload = edge_a.as_ref().map(|e| e.payload().read_arc());
    let expected = 2;
    let actual = edge_a_payload
      .as_ref()
      .and_then(|p| p.mutations.get("nuc"))
      .map_or(0, |v| v.len());
    assert_eq!(expected, actual);
    let expected = 1;
    let actual = edge_a_payload
      .as_ref()
      .and_then(|p| p.mutations.get("S"))
      .map_or(0, |v| v.len());
    assert_eq!(expected, actual);

    let output = auspice_write_str::<AuspiceWriter, _, _, _>(&graph)?;
    let graph2 = auspice_read_str::<AuspiceReader, _, _, _>(&output)?;

    assert!(graph2.data().read_arc().has_mutations);
    let expected = 3;
    let actual = graph2.num_nodes();
    assert_eq!(expected, actual);

    let edges2 = graph2.get_edges();
    let edge_a2 = edges2.iter().find(|e| {
      let edge = e.read_arc();
      let target_node = graph2.get_node(edge.target());
      target_node.is_some_and(|n| n.read_arc().payload().read_arc().name.as_deref() == Some("A"))
    });
    let edge_a2 = edge_a2.map(|e| e.read_arc());
    let edge_a2_payload = edge_a2.as_ref().map(|e| e.payload().read_arc());
    let expected = 2;
    let actual = edge_a2_payload
      .as_ref()
      .and_then(|p| p.mutations.get("nuc"))
      .map_or(0, |v| v.len());
    assert_eq!(expected, actual);
    let expected = 1;
    let actual = edge_a2_payload
      .as_ref()
      .and_then(|p| p.mutations.get("S"))
      .map_or(0, |v| v.len());
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_nwk_auspice_nwk_topology_preserved() -> Result<(), Report> {
    let nwk_input = "(A:1,B:2)root;";
    let graph = nwk_read_str::<ConverterNode, ConverterEdge, ConverterData>(nwk_input)?;

    let edges = graph.get_edges();
    let edge_a = edges.iter().find(|e| {
      let edge = e.read_arc();
      let target_node = graph.get_node(edge.target());
      target_node.is_some_and(|n| n.read_arc().payload().read_arc().name.as_deref() == Some("A"))
    });
    assert!(edge_a.is_none_or(|e| e.read_arc().payload().read_arc().mutations.is_empty()));

    let auspice_str = auspice_write_str::<AuspiceWriter, _, _, _>(&graph)?;
    let graph2 = auspice_read_str::<AuspiceReader, _, _, _>(&auspice_str)?;
    let nwk_output = nwk_write_str(&graph2, &NwkWriteOptions::default())?;

    assert_eq!(nwk_input, nwk_output);
    Ok(())
  }
}
