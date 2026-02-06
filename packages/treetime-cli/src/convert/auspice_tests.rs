use crate::convert::auspice::{AuspiceReader, AuspiceWriter};
use crate::convert::convert::{ConverterData, ConverterEdge, ConverterNode};
use crate::convert::mutation::PartitionedMutations;
use eyre::Report;
use indoc::indoc;
use pretty_assertions::assert_eq;
use treetime::graph::graph::Graph;
use treetime::io::auspice::{auspice_read_str, auspice_write_str};
use treetime::io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};
use treetime::make_internal_report;

fn get_mutations_for_node(
  graph: &Graph<ConverterNode, ConverterEdge, ConverterData>,
  node_name: &str,
) -> Result<PartitionedMutations, Report> {
  let node_key = graph
    .find_node(|n| n.name.as_deref() == Some(node_name))
    .ok_or_else(|| make_internal_report!("Node '{node_name}' not found"))?;
  let node = graph
    .get_node(node_key)
    .ok_or_else(|| make_internal_report!("Node key {node_key} not found"))?;
  let (_, edge) = graph.exactly_one_parent_of(&node.read_arc())?;
  Ok(edge.read_arc().payload().read_arc().mutations.clone())
}

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

  let mutations = get_mutations_for_node(&graph, "A")?;
  assert_eq!(2, mutations["nuc"].len());
  assert_eq!(1, mutations["S"].len());

  let output = auspice_write_str::<AuspiceWriter, _, _, _>(&graph)?;
  let graph = auspice_read_str::<AuspiceReader, _, _, _>(&output)?;

  assert!(graph.data().read_arc().has_mutations);
  assert_eq!(3, graph.num_nodes());

  let mutations = get_mutations_for_node(&graph, "A")?;
  assert_eq!(2, mutations["nuc"].len());
  assert_eq!(1, mutations["S"].len());

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
