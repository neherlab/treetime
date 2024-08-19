use crate::alphabet::alphabet::Alphabet;
use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named};
use crate::io::nwk::{EdgeFromNwk, NodeFromNwk};
use crate::port::composition::Composition;
use crate::port::mutation::{InDel, Mut};
use crate::port::seq_partitions::SeqPartition;
use crate::seq::find_char_ranges::find_letter_ranges;
use crate::seq::range_union::range_union;
use eyre::Report;
use maplit::btreemap;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type SparseGraph<'g> = Graph<SparseNode, SparseEdge, SparseMeta<'g>>;

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseNode {
  pub name: Option<String>,
  pub sparse_partitions: Vec<SparseSeqNode>,
}

impl NodeFromNwk for SparseNode {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|s| s.as_ref().to_owned()),
      ..SparseNode::default()
    })
  }
}

impl GraphNode for SparseNode {}

impl Named for SparseNode {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseSeqInfo {
  pub unknown: Vec<(usize, usize)>,
  pub gaps: Vec<(usize, usize)>,
  pub non_char: Vec<(usize, usize)>, // any position that does not evolve according to the substitution model, i.e. gap or N
  pub composition: Composition,      // count of all characters in the region that is not `non_char`
  pub sequence: Vec<char>,
  pub fitch: SparseSeqDis,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseSeqNode {
  pub seq: SparseSeqInfo,
  pub profile: SparseSeqDis,
  pub msg_to_parents: SparseSeqDis, // there might be multiple parents, but all parents only see info from children
  pub msgs_to_children: BTreeMap<String, SparseSeqDis>,
  pub msgs_from_children: BTreeMap<String, SparseSeqDis>,
}

impl SparseSeqNode {
  pub fn new(seq: &[char], alphabet: &Alphabet) -> Result<Self, Report> {
    // FIXME: the original code used `alphabet_gapN`:
    //
    // alphabet_gapN = ''.join(gtr.alphabet)+'-N'
    //
    // def seq_info_from_array(seq: np.array, profile: Dict[str,np.array], alphabet_gapN: str) -> SparseSeqInfo

    let variable = seq
      .iter()
      .enumerate()
      .filter(|(_, &c)| alphabet.is_ambiguous(c))
      .map(|(pos, &c)| {
        (
          pos,
          VarPos {
            dis: alphabet.get_profile(c).clone(),
            state: None,
          },
        )
      })
      .collect();

    let seq_dis = SparseSeqDis {
      variable,
      ..SparseSeqDis::default()
    };

    let unknown = find_letter_ranges(seq, alphabet.unknown());
    let gaps = find_letter_ranges(seq, alphabet.gap());
    let non_char = range_union(&[unknown.clone(), gaps.clone()]); // TODO(perf): avoid cloning

    Ok(Self {
      seq: SparseSeqInfo {
        unknown,
        gaps,
        non_char,
        composition: Composition::default(),
        sequence: seq.to_owned(), // TODO(perf): try to avoid cloning
        fitch: seq_dis,
      },
      profile: SparseSeqDis::default(),
      msg_to_parents: SparseSeqDis::default(),
      msgs_to_children: btreemap! {},
      msgs_from_children: btreemap! {},
    })
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseEdge {
  pub sparse_partitions: Vec<SparseSeqEdge>,
  pub branch_length: Option<f64>,
}

impl GraphEdge for SparseEdge {}

impl Weighted for SparseEdge {
  fn weight(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_weight(&mut self, weight: Option<f64>) {
    self.branch_length = weight;
  }
}

impl EdgeFromNwk for SparseEdge {
  fn from_nwk(_: Option<f64>) -> Result<Self, Report> {
    Ok(Self::default())
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseSeqEdge {
  pub muts: Vec<Mut>,
  pub indels: Vec<InDel>,
  pub transmission: Option<Vec<(usize, usize)>>,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseMeta<'g> {
  #[serde(skip)]
  pub sparse_partitions: Vec<SeqPartition<'g>>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct VarPos {
  pub dis: Array1<f64>, // array of floats of size 'alphabet'
  pub state: Option<char>,
}

impl VarPos {
  pub fn new(dis: Array1<f64>, state: Option<char>) -> Self {
    Self { dis, state }
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Deletion {
  pub deleted: usize, // number of times deletion is observed
  pub ins: usize,     // or not
  pub alt: Vec<char>,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseSeqDis {
  /// probability vector for each variable position collecting information from children
  pub variable: BTreeMap<usize, VarPos>,

  pub variable_indel: BTreeMap<(usize, usize), Deletion>,

  /// probability vector for the state of fixed positions based on information from children
  pub fixed: BTreeMap<char, Array1<f64>>,

  pub fixed_counts: Composition,

  /// Total log likelihood
  pub log_lh: f64,
}
