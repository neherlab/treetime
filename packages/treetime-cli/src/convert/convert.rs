use crate::convert::args::{Args, TreeFormat};
use crate::convert::auspice::{AuspiceReader, AuspiceWriter};
use crate::convert::usher::{UsherReader, UsherWriter};
use eyre::Report;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;
use treetime::graph::edge::{GraphEdge, Weighted};
use treetime::graph::graph::Graph;
use treetime::graph::node::{GraphNode, Named};
use treetime::io::auspice::{AuspiceTreeMeta, auspice_read_file, auspice_write_file};
use treetime::io::json::{JsonPretty, json_read_file, json_write_file};
use treetime::io::nex::{NexWriteOptions, nex_write_file};
use treetime::io::nwk::{NwkWriteOptions, nwk_read_file, nwk_write_file};
use treetime::io::phyloxml::{
  PhyloxmlJsonOptions, phyloxml_json_read_file, phyloxml_json_write_file, phyloxml_read_file, phyloxml_write_file,
};
use treetime::io::usher_mat::{
  UsherMatJsonOptions, usher_mat_json_read_file, usher_mat_json_write_file, usher_mat_pb_read_file,
  usher_mat_pb_write_file,
};

pub type ConverterGraph = Graph<ConverterNode, ConverterEdge, ConverterData>;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConverterNode {
  pub name: Option<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConverterEdge {
  pub weight: Option<f64>,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct ConverterData {
  pub rooted: bool,
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
