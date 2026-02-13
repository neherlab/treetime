use crate::convert::auspice::AuspiceReader;
use crate::convert::convert::{ConverterData, ConverterEdge, ConverterNode};
use crate::convert::usher::UsherWriter;
use eyre::Report;
use indoc::indoc;
use pretty_assertions::assert_eq;
use std::collections::BTreeSet;
use treetime_io::auspice::auspice_read_str;
use treetime_io::json::json_read_str;
use treetime_io::usher_mat::{UsherMatJsonOptions, UsherTree, usher_mat_json_write_str};

#[test]
fn test_usher_write_mutations() -> Result<(), Report> {
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
                  "nuc": ["A100T", "C200G"]
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

  let graph = auspice_read_str::<AuspiceReader, ConverterNode, ConverterEdge, ConverterData>(auspice_input)?;

  assert!(graph.data().read_arc().has_mutations);

  let usher_output = usher_mat_json_write_str::<UsherWriter, _, _, _>(&graph, &UsherMatJsonOptions::default())?;
  let usher_tree: UsherTree = json_read_str(&usher_output)?;

  let expected = BTreeSet::from([100, 200, 300]);
  let actual = usher_tree.get_all_positions();
  assert_eq!(expected, actual);

  Ok(())
}

#[test]
fn test_auspice_usher_auspice_nuc_mutations() -> Result<(), Report> {
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
                  "nuc": ["A100T", "C200G"],
                  "S": ["D614G"]
                }
              },
              "node_attrs": {
                "div": 2.0
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
  let expected = BTreeSet::from([100, 200]);
  assert_eq!(
    expected, positions,
    "Only nuc partition mutations should appear in UShER output"
  );

  Ok(())
}
