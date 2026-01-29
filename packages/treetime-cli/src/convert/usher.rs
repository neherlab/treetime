use crate::convert::convert::{ConverterData, ConverterEdge, ConverterNode};
use crate::convert::mutation::{Mutation, MutationList, PartitionedMutations};
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use treetime::graph::graph::Graph;
use treetime::io::usher_mat::{
  UsherGraphContext, UsherMetadata, UsherMutation, UsherMutationList, UsherRead, UsherTree, UsherTreeContext,
  UsherTreeNode, UsherWrite,
};
use treetime::make_internal_report;

const NUC_PARTITION: &str = "nuc";
const NUC_CHARS: [u8; 4] = [b'A', b'C', b'G', b'T'];

fn nuc_int_to_char(nuc: i32) -> Option<u8> {
  NUC_CHARS.get(nuc as usize).copied()
}

fn nuc_char_to_int(nuc: u8) -> Option<i32> {
  NUC_CHARS.iter().position(|&c| c == nuc).map(|i| i as i32)
}

fn parse_usher_mutations(usher_mutations: &[UsherMutation]) -> Result<PartitionedMutations, Report> {
  if usher_mutations.is_empty() {
    return Ok(btreemap! {});
  }

  let mutations: MutationList = usher_mutations
    .iter()
    .filter_map(|m| {
      let reference = nuc_int_to_char(m.par_nuc)?;
      let alternative = m.mut_nuc.first().and_then(|&n| nuc_int_to_char(n))?;
      let position = m.position as usize;
      Some(Mutation::new(reference, position, alternative))
    })
    .collect_vec();

  if mutations.is_empty() {
    Ok(btreemap! {})
  } else {
    Ok(btreemap! { NUC_PARTITION.to_owned() => mutations })
  }
}

fn format_usher_mutations(mutations: &PartitionedMutations) -> UsherMutationList {
  let nuc_mutations = mutations.get(NUC_PARTITION);

  let mutation = nuc_mutations.map_or(vec![], |muts| {
    muts
      .iter()
      .filter_map(|m| {
        let par_nuc = nuc_char_to_int(m.reference)?;
        let mut_nuc = nuc_char_to_int(m.alternative)?;
        Some(UsherMutation {
          position: m.position as i32,
          ref_nuc: par_nuc,
          par_nuc,
          mut_nuc: vec![mut_nuc],
          chromosome: String::new(),
        })
      })
      .collect_vec()
  });

  UsherMutationList { mutation }
}

pub struct UsherWriter;

impl UsherWrite<ConverterNode, ConverterEdge, ConverterData> for UsherWriter {
  fn new(_graph: &Graph<ConverterNode, ConverterEdge, ConverterData>) -> Result<Self, Report> {
    Ok(Self {})
  }

  fn usher_node_from_graph_components(
    &mut self,
    context: &UsherGraphContext<ConverterNode, ConverterEdge, ConverterData>,
  ) -> Result<(UsherTreeNode, UsherMutationList, UsherMetadata), Report> {
    let &UsherGraphContext { node, edge, .. } = context;

    let node_name = node
      .name
      .as_ref()
      .ok_or_else(|| make_internal_report!("Encountered node with empty name"))?
      .to_owned();

    let usher_node = UsherTreeNode {
      node_name,
      condensed_leaves: vec![],
    };

    let usher_mutations = edge.map_or(UsherMutationList { mutation: vec![] }, |e| {
      format_usher_mutations(&e.mutations)
    });

    let usher_meta = UsherMetadata {
      clade_annotations: vec![],
    };

    Ok((usher_node, usher_mutations, usher_meta))
  }
}

pub struct UsherReader {
  has_mutations: bool,
}

impl UsherRead<ConverterNode, ConverterEdge, ConverterData> for UsherReader {
  fn new(tree: &UsherTree) -> Result<Self, Report> {
    let has_mutations = tree.node_mutations.iter().any(|ml| !ml.mutation.is_empty());
    Ok(Self { has_mutations })
  }

  fn usher_data_to_graph_data(&mut self, _tree: &UsherTree) -> Result<ConverterData, Report> {
    Ok(ConverterData {
      has_mutations: self.has_mutations,
      ..Default::default()
    })
  }

  fn usher_node_to_graph_components(
    &mut self,
    context: &UsherTreeContext,
  ) -> Result<(ConverterNode, ConverterEdge), Report> {
    let UsherTreeContext { node, .. } = context;

    let mutations = parse_usher_mutations(&node.mutations)?;

    let node = ConverterNode {
      name: node.name.as_ref().map(ToOwned::to_owned),
    };

    let edge = ConverterEdge {
      weight: None,
      mutations,
    };

    Ok((node, edge))
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::convert::auspice::AuspiceReader;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeSet;
  use treetime::io::auspice::auspice_read_str;
  use treetime::io::usher_mat::{UsherMatJsonOptions, UsherTree, usher_mat_json_write_str};
  use treetime_io::json::json_read_str;

  fn get_all_positions(tree: &UsherTree) -> BTreeSet<i32> {
    tree
      .node_mutations
      .iter()
      .flat_map(|ml| ml.mutation.iter().map(|m| m.position))
      .collect()
  }

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
    let actual = get_all_positions(&usher_tree);
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

    let positions = get_all_positions(&usher_tree);
    let expected = BTreeSet::from([100, 200]);
    assert_eq!(
      expected, positions,
      "Only nuc partition mutations should appear in UShER output"
    );

    Ok(())
  }
}
