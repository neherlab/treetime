use crate::convert::args::{Args, TreeFormat};
use crate::convert::auspice::{AuspiceReader, AuspiceWriter};
use crate::convert::mutation::PartitionedMutations;
use crate::convert::usher::{UsherReader, UsherWriter};
use eyre::Report;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;
use treetime::graph::edge::{GraphEdge, Weighted};
use treetime::graph::graph::Graph;
use treetime::graph::node::{GraphNode, Named};
use treetime::io::auspice::{AuspiceTreeMeta, auspice_read_file, auspice_write_file};
use treetime::io::nex::{NexWriteOptions, nex_write_file};
use treetime::io::nwk::{NwkWriteOptions, nwk_read_file, nwk_write_file};
use treetime::io::phyloxml::{
  PhyloxmlJsonOptions, phyloxml_json_read_file, phyloxml_json_write_file, phyloxml_read_file, phyloxml_write_file,
};
use treetime::io::usher_mat::{
  UsherMatJsonOptions, usher_mat_json_read_file, usher_mat_json_write_file, usher_mat_pb_read_file,
  usher_mat_pb_write_file,
};
use treetime_io::json::{JsonPretty, json_read_file, json_write_file};

pub type ConverterGraph = Graph<ConverterNode, ConverterEdge, ConverterData>;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConverterNode {
  pub name: Option<String>,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct ConverterEdge {
  pub weight: Option<f64>,
  #[serde(default, skip_serializing_if = "PartitionedMutations::is_empty")]
  pub mutations: PartitionedMutations,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct ConverterData {
  pub rooted: bool,
  pub has_mutations: bool,
  pub version: Option<String>,
  pub meta: AuspiceTreeMeta,
  pub root_sequence: Option<BTreeMap<String, String>>,
  pub other: Value,
}

impl GraphNode for ConverterNode {}

impl Named for ConverterNode {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

impl GraphEdge for ConverterEdge {}

impl Weighted for ConverterEdge {
  fn weight(&self) -> Option<f64> {
    self.weight
  }
  fn set_weight(&mut self, weight: Option<f64>) {
    self.weight = weight;
  }
}

pub fn converter_write_file(args: &Args, output_format: TreeFormat, graph: &ConverterGraph) -> Result<(), Report> {
  match output_format {
    TreeFormat::Auspice => auspice_write_file::<AuspiceWriter, _, _, _>(&args.output, graph),
    TreeFormat::Newick => nwk_write_file(&args.output, graph, &NwkWriteOptions::default()),
    TreeFormat::Nexus => nex_write_file(&args.output, graph, &NexWriteOptions::default()),
    TreeFormat::PhyloGraph => json_write_file(&args.output, graph, JsonPretty(true)),
    TreeFormat::MatJson => {
      usher_mat_json_write_file::<UsherWriter, _, _, _>(&args.output, graph, &UsherMatJsonOptions::default())
    },
    TreeFormat::MatPb => usher_mat_pb_write_file::<UsherWriter, _, _, _>(&args.output, graph),
    TreeFormat::Phyloxml => phyloxml_write_file(&args.output, graph),
    TreeFormat::PhyloxmlJson => phyloxml_json_write_file(&args.output, graph, &PhyloxmlJsonOptions::default()),
  }
}

pub fn converter_read_file(args: &Args, input_format: TreeFormat) -> Result<ConverterGraph, Report> {
  match input_format {
    TreeFormat::Auspice => auspice_read_file::<AuspiceReader, _, _, _>(&args.input),
    TreeFormat::Newick => nwk_read_file(&args.input),
    TreeFormat::Nexus => unimplemented!("Reading Nexus files is not yet implemented"),
    TreeFormat::PhyloGraph => json_read_file(&args.input),
    TreeFormat::MatJson => usher_mat_json_read_file::<UsherReader, _, _, _>(&args.input),
    TreeFormat::MatPb => usher_mat_pb_read_file::<UsherReader, _, _, _>(&args.input),
    TreeFormat::Phyloxml => phyloxml_read_file(&args.input),
    TreeFormat::PhyloxmlJson => phyloxml_json_read_file(&args.input),
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::convert::auspice::AuspiceReader;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeSet;
  use treetime::io::auspice::{auspice_read_str, auspice_write_str};
  use treetime::io::nwk::{nwk_read_str, nwk_write_str};
  use treetime::io::usher_mat::{UsherMatJsonOptions, UsherTree, usher_mat_json_write_str};
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
}
