use crate::convert::convert::{ConverterData, ConverterEdge, ConverterNode};
use eyre::Report;
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

    Ok(AuspiceTreeNode {
      name,
      branch_attrs: AuspiceTreeBranchAttrs::default(),
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

pub struct AuspiceReader;

impl AuspiceRead<ConverterNode, ConverterEdge, ConverterData> for AuspiceReader {
  fn new(_: &AuspiceTree) -> Result<Self, Report> {
    Ok(Self {})
  }

  fn auspice_data_to_graph_data(&mut self, tree: &AuspiceTree) -> Result<ConverterData, Report> {
    Ok(ConverterData {
      rooted: true,
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
    Ok((
      ConverterNode {
        name: Some(node.name.clone()),
      },
      ConverterEdge { weight: branch_length },
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
}
