use crate::alphabet::alphabet::Alphabet;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named};
use crate::port::mutation::{InDel, Mut};
use crate::port::seq_partitions::SeqPartition;
use crate::seq::find_char_ranges::find_letter_ranges;
use eyre::Report;
use maplit::btreemap;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type SparseSeqGraph<'a, 'g> = Graph<SparseSeqNode, SparseSeqEdge, SparseSeqMeta<'a, 'g>>;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct VarPos {
  dis: Array1<f64>, // array of floats of size 'alphabet'
  state: Option<char>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Deletion {
  deleted: usize, // number of times deletion is observed
  ins: usize,     // or not
  alt: String,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseSeqDis {
  /// probability vector for each variable position collecting information from children
  variable: BTreeMap<usize, VarPos>,

  variable_indel: BTreeMap<(usize, usize), Deletion>,

  /// probability vector for the state of fixed positions based on information from children
  fixed: BTreeMap<String, Array1<f64>>,

  fixed_counts: BTreeMap<String, usize>,

  /// Total log likelihood
  log_lh: f64,
}

// #[derive(Debug)]
// pub struct FitchVar {
//   /// Probability vector for each variable position collecting information from children
//   variable: BTreeMap<usize, VarPos>,
//   variable_indel: BTreeMap<(usize, usize), Deletion>,
// }

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseSeqInfo {
  unknown: Vec<(usize, usize)>,
  gaps: Vec<(usize, usize)>,
  non_char: Vec<(usize, usize)>, // any position that does not evolve according to the substitution model, i.e. gap or N
  composition: BTreeMap<String, usize>, // count of all characters in the region that is not `non_char`
  sequence: Option<Array1<f64>>,
  fitch: SparseSeqDis,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseSeqNode {
  pub name: Option<String>,
  pub seq: SparseSeqInfo,
  pub profile: SparseSeqDis,
  pub msg_to_parents: SparseSeqDis, // there might be multiple parents, but all parents only see info from children
  pub msgs_to_children: BTreeMap<String, SparseSeqDis>,
  pub msgs_from_children: BTreeMap<String, SparseSeqDis>,
}

impl GraphNode for SparseSeqNode {}

impl Named for SparseSeqNode {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

impl SparseSeqNode {
  pub fn new(name: String, seq: &[char], alphabet: &Alphabet) -> Result<Self, Report> {
    // FIXME: the original code used `alphabet_gapN`:
    //
    // alphabet_gapN = ''.join(gtr.alphabet)+'-N'
    //
    // def seq_info_from_array(seq: np.array, profile: Dict[str,np.array], alphabet_gapN: str) -> SparseSeqInfo

    let variable = seq
      .iter()
      .filter(|&c| alphabet.contains(*c) && !(alphabet.is_gap(*c) || alphabet.is_ambiguous(*c)))
      .enumerate()
      .map(|(pos, c)| {
        (
          pos,
          VarPos {
            dis: alphabet.get_profile(*c).clone(),
            state: None,
          },
        )
      })
      .collect();

    let seq_dis = SparseSeqDis {
      variable,
      ..SparseSeqDis::default()
    };

    let unknown = find_letter_ranges(seq, alphabet.ambiguous());

    let gaps = alphabet
      .gap()
      .map(|gap| find_letter_ranges(seq, gap))
      .unwrap_or_default();

    Ok(Self {
      name: Some(name),
      seq: SparseSeqInfo {
        unknown,
        gaps,
        non_char: vec![],
        composition: btreemap! {},
        sequence: None,
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
pub struct SparseSeqEdge {
  muts: Vec<Mut>,
  indels: Vec<InDel>,
  transmission: Option<Vec<(usize, usize)>>,
}

impl GraphEdge for SparseSeqEdge {}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseSeqMeta<'a, 'g> {
  #[serde(skip)]
  pub sparse_partitions: Vec<SeqPartition<'a, 'g>>,
}
