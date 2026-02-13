use crate::convert::auspice::{AuspiceReader, AuspiceWriter};
use crate::convert::convert::{ConverterData, ConverterEdge, ConverterNode};
use crate::convert::usher::UsherWriter;
use eyre::Report;
use indoc::indoc;
use pretty_assertions::assert_eq;
use std::collections::BTreeSet;
use treetime_io::auspice::{auspice_read_str, auspice_write_str};
use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};
use treetime_io::usher_mat::{UsherMatJsonOptions, UsherTree, usher_mat_json_write_str};
use treetime_io::json::json_read_str;

#[test]
fn test_newick_to_auspice_to_newick_topology() -> Result<(), Report> {
  let nwk_input = "(A:1,B:2)root;";
  let graph = nwk_read_str::<ConverterNode, ConverterEdge, ConverterData>(nwk_input)?;

  let expected = 3;
  let actual = graph.num_nodes();
  assert_eq!(expected, actual);

  let auspice_str = auspice_write_str::<AuspiceWriter, _, _, _>(&graph)?;
  let graph2 = auspice_read_str::<AuspiceReader, ConverterNode, ConverterEdge, ConverterData>(&auspice_str)?;

  let expected = 3;
  let actual = graph2.num_nodes();
  assert_eq!(expected, actual);

  let nwk_output = nwk_write_str(&graph2, &NwkWriteOptions::default())?;
  assert_eq!(nwk_input, nwk_output);

  Ok(())
}

#[test]
fn test_auspice_to_newick_to_auspice_topology() -> Result<(), Report> {
  let auspice_input = indoc!(
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
              "node_attrs": {
                "div": 1.0
              }
            },
            {
              "name": "B",
              "node_attrs": {
                "div": 2.0
              }
            }
          ]
        }
      }"#
  );

  let graph = auspice_read_str::<AuspiceReader, ConverterNode, ConverterEdge, ConverterData>(auspice_input)?;

  let expected = 3;
  let actual = graph.num_nodes();
  assert_eq!(expected, actual);

  let nwk_str = nwk_write_str(&graph, &NwkWriteOptions::default())?;
  let graph2 = nwk_read_str::<ConverterNode, ConverterEdge, ConverterData>(&nwk_str)?;

  let expected = 3;
  let actual = graph2.num_nodes();
  assert_eq!(expected, actual);

  Ok(())
}

#[test]
fn test_auspice_to_usher_mutation_transfer() -> Result<(), Report> {
  let auspice_input = indoc!(
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
                  "nuc": ["A100T"],
                  "S": ["D614G"]
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

  let graph = auspice_read_str::<AuspiceReader, ConverterNode, ConverterEdge, ConverterData>(auspice_input)?;
  let usher_output = usher_mat_json_write_str::<UsherWriter, _, _, _>(&graph, &UsherMatJsonOptions::default())?;
  let usher_tree: UsherTree = json_read_str(&usher_output)?;

  let positions = usher_tree.get_all_positions();
  let expected = BTreeSet::from([100]);
  assert_eq!(
    expected, positions,
    "Only nuc partition mutations should appear in UShER output"
  );

  Ok(())
}
