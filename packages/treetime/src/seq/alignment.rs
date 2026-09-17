use crate::make_error;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AlignmentRecord, Seq};

/// The merged tree-and-alignment input the reconstruction pipeline consumes: the graph, the per-node
/// input keyed by the graph's own node keys, and the per-edge input keyed by its edge keys.
#[derive(Debug)]
pub struct AncestralInput {
  pub graph: Graph,
  pub nodes: BTreeMap<GraphNodeKey, NodeSeqInput>,
  pub edges: BTreeMap<GraphEdgeKey, EdgeSeqInput>,
}

impl AncestralInput {
  /// The per-node name map, keyed by node id, as `assign_node_names`, the topology loops, and the
  /// output writers consume it.
  pub fn names(&self) -> BTreeMap<GraphNodeKey, Option<String>> {
    self.nodes.iter().map(|(key, node)| (*key, node.name.clone())).collect()
  }

  /// The per-edge branch-length map, keyed by edge id, as the reconstruction and output writers
  /// consume it.
  pub fn branch_lengths(&self) -> BTreeMap<GraphEdgeKey, Option<f64>> {
    self
      .edges
      .iter()
      .map(|(key, edge)| (*key, edge.branch_length))
      .collect()
  }
}

/// One node's reconstruction input: its name and its attached alignment sequence. `name` is `None`
/// where a node has no name; `seq` is `None` for internal nodes and for leaves absent from the
/// alignment.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NodeSeqInput {
  pub name: Option<String>,
  pub seq: Option<Seq>,
}

/// One edge's reconstruction input: the raw input-tree branch length, `None` where the edge carried
/// no `:length`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct EdgeSeqInput {
  pub branch_length: Option<f64>,
}

/// Build the per-node reconstruction input map from node names and alignment records, matching each
/// leaf to its record by name (first record wins on duplicate names).
///
/// Leaves with no matching record get `seq = None`; callers that require a sequence for every leaf
/// complete the alignment first. Internal nodes always get `seq = None`. Extra records that match no
/// leaf are ignored.
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

/// The common length of the attached leaf sequences in a node-input map, or `0` when none carry a
/// sequence. Errors when leaves disagree on length, listing each length and its leaf names.
pub fn get_common_length_of_node_inputs(node_inputs: &BTreeMap<GraphNodeKey, NodeSeqInput>) -> Result<usize, Report> {
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
