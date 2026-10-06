use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::mask::create_mask;
use crate::error::input_error;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pair_by_name::pair_by_name;
use treetime_primitives::{AlignmentRecord, Seq};

#[derive(Debug)]
pub struct AncestralInput {
  pub graph: Graph,
  pub nodes: BTreeMap<GraphNodeKey, NodeSeqInput>,
  pub edges: BTreeMap<GraphEdgeKey, EdgeSeqInput>,
  pub alphabet: Alphabet,
  pub mask: Vec<bool>,
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

#[derive(Clone, Debug, Default, PartialEq)]
pub struct EdgeSeqInput {
  pub branch_length: Option<f64>,
}

pub fn pair_leaf_sequences(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  records: Vec<AlignmentRecord>,
) -> LeafPairing {
  let pairing = pair_by_name(
    graph.get_leaves().map(|leaf| leaf.key()),
    names,
    records.into_iter().map(|record| (record.name, record.seq)),
  );
  LeafPairing {
    sequences: LeafSequences::new(names, pairing.by_node, pairing.unmatched),
    duplicate_names: pairing.duplicate_entry_names,
  }
}

pub struct LeafPairing {
  pub sequences: LeafSequences,
  pub duplicate_names: Vec<String>,
}

#[derive(Clone, Debug)]
pub struct LeafSequences {
  pub nodes: BTreeMap<GraphNodeKey, NodeSeqInput>,
  pub unmatched: Vec<(String, Seq)>,
}

impl LeafSequences {
  pub fn new(
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    mut seqs: BTreeMap<GraphNodeKey, Seq>,
    unmatched: Vec<(String, Seq)>,
  ) -> Self {
    let nodes = names
      .iter()
      .map(|(key, name)| {
        let input = NodeSeqInput {
          name: name.clone(),
          seq: seqs.remove(key),
        };
        (*key, input)
      })
      .collect();
    Self { nodes, unmatched }
  }

  pub fn common_length(&self) -> Result<usize, Report> {
    common_length(self.kept())
  }

  pub fn mask(&self, alignment_length: usize, alphabet: &Alphabet) -> Vec<bool> {
    create_mask(self.kept().map(|(_, seq)| seq), alignment_length, alphabet)
  }

  fn kept(&self) -> impl Iterator<Item = (&str, &Seq)> {
    let paired = self
      .nodes
      .values()
      .filter_map(|node| Some((node.name.as_deref().unwrap_or(""), node.seq.as_ref()?)));
    let unmatched = self.unmatched.iter().map(|(name, seq)| (name.as_str(), seq));
    paired.chain(unmatched)
  }
}

pub(crate) fn get_common_length_of_node_inputs(
  node_inputs: &BTreeMap<GraphNodeKey, NodeSeqInput>,
) -> Result<usize, Report> {
  common_length(
    node_inputs
      .values()
      .filter_map(|node| Some((node.name.as_deref().unwrap_or(""), node.seq.as_ref()?))),
  )
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NodeSeqInput {
  pub name: Option<String>,
  pub seq: Option<Seq>,
}

pub fn get_common_length(aln: &[AlignmentRecord]) -> Result<usize, Report> {
  common_length(aln.iter().map(|record| (record.name.as_str(), &record.seq)))
}

fn common_length<'a>(seqs: impl IntoIterator<Item = (&'a str, &'a Seq)>) -> Result<usize, Report> {
  let lengths = seqs
    .into_iter()
    .into_group_map_by(|(_, seq)| seq.len())
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
          let names = entries.iter().map(|(name, _)| format!("    \"{name}\"")).join("\n");
          format!("Length {length}:\n{names}")
        })
        .join("\n\n");

      Err(input_error(format!(
        "Sequences are expected to all have the same length, but found the following lengths:\n\n{message}"
      )))
    },
  }
  .wrap_err("When calculating length of sequences")
}
