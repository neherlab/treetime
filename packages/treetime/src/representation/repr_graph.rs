use crate::graph::edge::{GraphEdge, NumMuts, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{Described, GraphNode, Named};
use crate::io::graphviz::{EdgeToGraphViz, NodeToGraphviz};
use crate::io::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions, format_weight};
use crate::o;
use crate::representation::edge_partition::EdgePartition;
use crate::representation::node_partition::NodePartition;
use eyre::Report;
use std::collections::BTreeMap;

pub type ReprGraph = Graph<ReprNode, ReprEdge, ()>;

#[derive(Default, Debug)]
pub struct ReprNode {
  name: Option<String>,
  desc: Option<String>,
  partitions: Vec<Box<dyn NodePartition>>,
}

impl ReprNode {
  pub fn new(name: Option<String>, desc: Option<String>, partitions: Vec<Box<dyn NodePartition>>) -> Self {
    Self { name, desc, partitions }
  }

  pub fn partitions(&self) -> &[Box<dyn NodePartition>] {
    &self.partitions
  }

  pub fn set_partitions(&mut self, partitions: Vec<Box<dyn NodePartition>>) {
    self.partitions = partitions;
  }

  pub fn partition_at(&self, index: usize) -> &Box<dyn NodePartition> {
    &self.partitions[index]
  }

  pub fn set_partition_at(&mut self, index: usize, partition: Box<dyn NodePartition>) {
    self.partitions[index] = partition;
  }
}

impl NodeFromNwk for ReprNode {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|s| s.as_ref().to_owned()),
      ..ReprNode::default()
    })
  }
}

impl NodeToNwk for ReprNode {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    let mutations: String = "".to_owned(); // TODO: fill mutations
    BTreeMap::from([(o!("mutations"), mutations)])
  }
}

impl GraphNode for ReprNode {}

impl Named for ReprNode {
  fn get_name_maybe(&self) -> Option<&str> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<String>) {
    self.name = name;
  }
}

impl Described for ReprNode {
  fn get_desc(&self) -> Option<&str> {
    self.desc.as_deref()
  }

  fn set_desc(&mut self, desc: Option<String>) {
    self.desc = desc;
  }
}

impl NodeToGraphviz for ReprNode {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

#[derive(Default, Debug)]
pub struct ReprEdge {
  branch_length: Option<f64>,
  partitions: Vec<Box<dyn EdgePartition>>,
}

impl ReprEdge {
  pub fn new(branch_length: Option<f64>, partitions: Vec<Box<dyn EdgePartition>>) -> Self {
    Self {
      branch_length,
      partitions,
    }
  }

  pub fn partitions(&self) -> &[Box<dyn EdgePartition>] {
    &self.partitions
  }

  pub fn set_partitions(&mut self, partitions: Vec<Box<dyn EdgePartition>>) {
    self.partitions = partitions;
  }

  pub fn partition_at(&self, index: usize) -> &Box<dyn EdgePartition> {
    &self.partitions[index]
  }

  pub fn set_partition_at(&mut self, index: usize, partition: Box<dyn EdgePartition>) {
    self.partitions[index] = partition;
  }
}

impl GraphEdge for ReprEdge {}

impl Weighted for ReprEdge {
  fn weight(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_weight(&mut self, weight: Option<f64>) {
    self.branch_length = weight;
  }
}

impl NumMuts for ReprEdge {
  fn num_muts(&self) -> Option<usize> {
    if self.partitions.is_empty() {
      None // Unknown number of mutations when no partitions are present
    } else {
      Some(self.partitions.iter().map(|p| p.subs.len()).sum())
    }
  }
}

impl EdgeFromNwk for ReprEdge {
  fn from_nwk(branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self {
      branch_length,
      ..Self::default()
    })
  }
}

impl EdgeToNwk for ReprEdge {
  fn nwk_weight(&self) -> Option<f64> {
    self.weight()
  }
}

impl EdgeToGraphViz for ReprEdge {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .weight()
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.weight()
  }
}
