use crate::alphabet::alphabet::Alphabet;
use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named};
use crate::io::nwk::{EdgeFromNwk, NodeFromNwk};
use crate::port::mutation::InDel;
use crate::port::seq_partitions::SeqPartition;
use crate::seq::find_char_ranges::find_letter_ranges;
use eyre::Report;
use maplit::btreemap;
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type DenseGraph = Graph<DenseNode, DenseEdge, ()>;

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseNode {
  pub name: Option<String>,
  pub dense_partitions: Vec<DenseSeqNode>,
}

impl NodeFromNwk for DenseNode {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|s| s.as_ref().to_owned()),
      ..DenseNode::default()
    })
  }
}

impl GraphNode for DenseNode {}

impl Named for DenseNode {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseSeqInfo {
  pub gaps: Vec<(usize, usize)>,
  pub sequence: Vec<char>,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseSeqNode {
  pub seq: DenseSeqInfo,
  pub profile: DenseSeqDis,
  pub msg_to_parents: DenseSeqDis, // there might be multiple parents, but all parents only see info from children
  pub msgs_to_children: BTreeMap<String, DenseSeqDis>,
  pub msgs_from_children: BTreeMap<String, DenseSeqDis>,
}

impl DenseSeqNode {
  pub fn new(seq: &[char], alphabet: &Alphabet) -> Result<Self, Report> {
    let gaps = find_letter_ranges(seq, alphabet.gap());

    Ok(Self {
      seq: DenseSeqInfo {
        gaps,
        sequence: seq.to_owned(), // TODO(perf): try to avoid cloning
      },
      profile: DenseSeqDis::default(),
      msg_to_parents: DenseSeqDis::new(alphabet.seq2prof(seq)?),
      msgs_to_children: btreemap! {},
      msgs_from_children: btreemap! {},
    })
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseEdge {
  pub dense_partitions: Vec<DenseSeqEdge>,
  pub branch_length: Option<f64>,
}

impl GraphEdge for DenseEdge {}

impl Weighted for DenseEdge {
  fn weight(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_weight(&mut self, weight: Option<f64>) {
    self.branch_length = weight;
  }
}

impl EdgeFromNwk for DenseEdge {
  fn from_nwk(branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self {
      branch_length,
      ..Self::default()
    })
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseSeqEdge {
  pub indels: Vec<InDel>,
  pub transmission: Option<Vec<(usize, usize)>>,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseMeta {
  #[serde(skip)]
  pub dense_partitions: Vec<SeqPartition>,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseSeqDis {
  pub dis: Array2<f64>,

  /// Total log likelihood
  pub log_lh: f64,
}

impl DenseSeqDis {
  pub fn new(dis: Array2<f64>) -> Self {
    Self { dis, log_lh: 0.0 }
  }
}
