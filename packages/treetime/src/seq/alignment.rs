use crate::make_error;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AlignmentRecord, Seq};

#[derive(Debug)]
pub struct AncestralInput {
  pub graph: Graph,
  pub nodes: BTreeMap<GraphNodeKey, NodeSeqInput>,
  pub edges: BTreeMap<GraphEdgeKey, EdgeSeqInput>,
}

impl AncestralInput {
  pub fn names(&self) -> BTreeMap<GraphNodeKey, Option<String>> {
    self.nodes.iter().map(|(key, node)| (*key, node.name.clone())).collect()
  }

  pub fn branch_lengths(&self) -> BTreeMap<GraphEdgeKey, Option<f64>> {
    self
      .edges
      .iter()
      .map(|(key, edge)| (*key, edge.branch_length))
      .collect()
  }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NodeSeqInput {
  pub name: Option<String>,
  pub seq: Option<Seq>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct EdgeSeqInput {
  pub branch_length: Option<f64>,
}

pub fn node_seq_inputs(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  aln: Vec<AlignmentRecord>,
) -> BTreeMap<GraphNodeKey, NodeSeqInput> {
  let mut records_by_name: BTreeMap<String, AlignmentRecord> = BTreeMap::new();
  for record in aln {
    records_by_name.entry(record.name.clone()).or_insert(record);
  }

  let leaf_keys: BTreeSet<GraphNodeKey> = graph.get_leaves().map(|leaf| leaf.key()).collect();

  names
    .iter()
    .map(|(key, name)| {
      let seq = if leaf_keys.contains(key) {
        name
          .as_deref()
          .and_then(|name| records_by_name.remove(name))
          .map(|record| record.seq)
      } else {
        None
      };
      (
        *key,
        NodeSeqInput {
          name: name.clone(),
          seq,
        },
      )
    })
    .collect()
}

pub(crate) fn get_common_length_of_node_inputs(
  node_inputs: &BTreeMap<GraphNodeKey, NodeSeqInput>,
) -> Result<usize, Report> {
  let lengths = node_inputs
    .values()
    .filter_map(|node| node.seq.as_ref().map(|seq| (seq.len(), node)))
    .into_group_map_by(|(length, _)| *length)
    .into_iter()
    .collect_vec();

  match lengths[..] {
    [] => Ok(0),
    [(length, _)] => Ok(length),
    _ => {
      let message = lengths
        .into_iter()
        .sorted_by_key(|(length, _)| *length)
        .map(|(length, entries)| {
          let names = entries
            .iter()
            .map(|(_, node)| format!("    \"{}\"", node.name.as_deref().unwrap_or("")))
            .join("\n");
          format!("Length {length}:\n{names}")
        })
        .join("\n\n");

      make_error!("Sequences are expected to all have the same length, but found the following lengths:\n\n{message}")
    },
  }
  .wrap_err("When calculating length of sequences")
}

pub fn get_common_length(aln: &[AlignmentRecord]) -> Result<usize, Report> {
  let lengths = aln
    .iter()
    .into_group_map_by(|aln| aln.seq.len())
    .into_iter()
    .collect_vec();

  match lengths[..] {
    [] => Ok(0),
    [(length, _)] => Ok(length),
    _ => {
      let message = lengths
        .into_iter()
        .sorted_by_key(|(length, _)| *length)
        .map(|(length, entries)| {
          let names = entries.iter().map(|aln| format!("    \"{}\"", aln.name)).join("\n");
          format!("Length {length}:\n{names}")
        })
        .join("\n\n");

      make_error!("Sequences are expected to all have the same length, but found the following lengths:\n\n{message}")
    },
  }
  .wrap_err("When calculating length of sequences")
}
