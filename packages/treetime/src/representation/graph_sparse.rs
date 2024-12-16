use crate::alphabet::alphabet::Alphabet;
use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named};
use crate::io::graphviz::{EdgeToGraphViz, NodeToGraphviz};
use crate::io::nwk::{format_weight, EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions};
use crate::o;
use crate::representation::seq::Seq;
use crate::representation::seq_char::AsciiChar;
use crate::representation::state_set::StateSet;
use crate::seq::composition::Composition;
use crate::seq::find_char_ranges::find_letter_ranges;
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use crate::utils::interval::range_union::range_union;
use eyre::Report;
use maplit::btreemap;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type SparseGraph = Graph<SparseNode, SparseEdge, ()>;

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseNode {
  pub name: Option<String>,
  pub desc: Option<String>,
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

impl NodeToNwk for SparseNode {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    let mutations: String = "".to_owned(); // TODO: fill mutations
    BTreeMap::from([(o!("mutations"), mutations)])
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

impl NodeToGraphviz for SparseNode {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseSeqInfo {
  pub unknown: Vec<(usize, usize)>,
  pub gaps: Vec<(usize, usize)>,
  pub non_char: Vec<(usize, usize)>, // any position that does not evolve according to the substitution model, i.e. gap or N
  pub composition: Composition,      // count of all characters in the region that is not `non_char`
  pub sequence: Seq,
  pub fitch: ParsimonySeqDis,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseSeqNode {
  pub seq: SparseSeqInfo,
  pub profile: SparseSeqDis,
}

impl SparseSeqNode {
  pub fn new(seq: &Seq, alphabet: &Alphabet) -> Result<Self, Report> {
    // FIXME: the original code used `alphabet_gapN`:
    //
    // alphabet_gapN = ''.join(gtr.alphabet)+'-N'
    //
    // def seq_info_from_array(seq: np.array, profile: Dict[str,np.array], alphabet_gapN: str) -> SparseSeqInfo

    let variable = seq
      .iter()
      .enumerate()
      .filter(|(_, &c)| alphabet.is_ambiguous(c))
      .map(|(pos, &c)| (pos, alphabet.char_to_set(c)))
      .collect();

    let seq_dis = ParsimonySeqDis {
      variable,
      variable_indel: btreemap! {},
      composition: Composition::new(alphabet.chars(), alphabet.gap()),
    };

    let unknown = find_letter_ranges(seq, alphabet.unknown());
    let gaps = find_letter_ranges(seq, alphabet.gap());
    let non_char = range_union(&[unknown.clone(), gaps.clone()]); // TODO(perf): avoid cloning

    Ok(Self {
      seq: SparseSeqInfo {
        unknown,
        gaps,
        non_char,
        composition: Composition::new(alphabet.chars(), alphabet.gap()),
        sequence: seq.to_owned(), // TODO(perf): try to avoid cloning
        fitch: seq_dis,
      },
      profile: SparseSeqDis {
        variable: btreemap! {},
        variable_indel: btreemap! {},
        fixed: btreemap! {},
        fixed_counts: Composition::new(alphabet.chars(), alphabet.gap()),
        log_lh: 0.0,
      },
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
  fn from_nwk(branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self {
      branch_length,
      ..Self::default()
    })
  }
}

impl EdgeToNwk for SparseEdge {
  fn nwk_weight(&self) -> Option<f64> {
    self.weight()
  }
}

impl EdgeToGraphViz for SparseEdge {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .weight()
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.weight()
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseSeqEdge {
  pub subs: Vec<Sub>,
  pub indels: Vec<InDel>,
  pub msg_to_parent: SparseSeqDis,
  pub msg_to_child: SparseSeqDis,
  pub msg_from_child: SparseSeqDis,
  pub transmission: Option<Vec<(usize, usize)>>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct VarPos {
  pub dis: Array1<f64>, // array of floats of size 'alphabet'
  pub state: AsciiChar,
}

impl VarPos {
  pub fn new(dis: Array1<f64>, state: AsciiChar) -> Self {
    Self { dis, state }
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Deletion {
  pub deleted: usize, // number of times deletion is observed
  pub present: usize, // or not
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseSeqDis {
  /// probability vector for each variable position collecting information from children
  pub variable: BTreeMap<usize, VarPos>,

  pub variable_indel: BTreeMap<(usize, usize), Deletion>,

  /// probability vector for the state of fixed positions based on information from children
  pub fixed: BTreeMap<AsciiChar, Array1<f64>>,

  pub fixed_counts: Composition,

  /// Total log likelihood
  pub log_lh: f64,
}

impl Default for SparseSeqDis {
  fn default() -> Self {
    Self {
      variable: btreemap! {},
      variable_indel: btreemap! {},
      fixed: btreemap! {},
      fixed_counts: Composition::new(std::iter::empty::<u8>(), b'-'),
      log_lh: 0.0,
    }
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ParsimonySeqDis {
  /// probability vector for each variable position collecting information from children
  pub variable: BTreeMap<usize, StateSet>,

  pub variable_indel: BTreeMap<(usize, usize), Deletion>,

  pub composition: Composition,
}
