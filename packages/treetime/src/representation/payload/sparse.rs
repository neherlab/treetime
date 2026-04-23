use crate::alphabet::alphabet::Alphabet;
use crate::seq::composition::Composition;
use crate::seq::find_char_ranges::find_letter_ranges;
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use eyre::Report;
use maplit::btreemap;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use treetime_primitives::AlphabetLike;
use treetime_primitives::{AsciiChar, Seq, StateSet, seq};
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseNodePartition {
  pub seq: SparseSeqInfo,
  pub profile: MarginalSparseSeqDistribution,
}

impl SparseNodePartition {
  /// Create an empty placeholder node partition for edge split operations.
  /// The sequence and profile will be computed during the marginal update pass.
  pub fn empty(alphabet: &Alphabet) -> Self {
    Self {
      seq: SparseSeqInfo {
        unknown: vec![],
        gaps: vec![],
        non_char: vec![],
        composition: Composition::new(alphabet.chars(), alphabet.gap()),
        sequence: seq![],
        fitch: FitchSeqDistribution {
          variable: btreemap! {},
          variable_indel: btreemap! {},
          composition: Composition::new(alphabet.chars(), alphabet.gap()),
          chosen_state: btreemap! {},
        },
      },
      profile: MarginalSparseSeqDistribution::default(),
    }
  }

  pub fn new(seq: &Seq, alphabet: &Alphabet) -> Result<Self, Report> {
    let variable = seq
      .iter()
      .enumerate()
      .filter(|&(_, c)| alphabet.is_ambiguous(*c))
      .map(|(pos, &c)| (pos, alphabet.char_to_set(c)))
      .collect();

    let seq_dis = FitchSeqDistribution {
      variable,
      variable_indel: btreemap! {},
      composition: Composition::new(alphabet.chars(), alphabet.gap()),
      chosen_state: btreemap! {},
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
      profile: MarginalSparseSeqDistribution {
        variable: btreemap! {},
        variable_indel: btreemap! {},
        fixed: btreemap! {},
        fixed_counts: Composition::new(alphabet.chars(), alphabet.gap()),
        log_lh: 0.0,
      },
    })
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseSeqInfo {
  pub unknown: Vec<(usize, usize)>,
  pub gaps: Vec<(usize, usize)>,
  pub non_char: Vec<(usize, usize)>, // any position that does not evolve according to the substitution model, i.e. gap or N
  pub composition: Composition,      // count of all characters in the region that is not `non_char`
  pub sequence: Seq,
  pub fitch: FitchSeqDistribution,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseEdgePartition {
  pub subs: Vec<Sub>,
  pub indels: Vec<InDel>,
  pub msg_to_parent: MarginalSparseSeqDistribution,
  pub msg_to_child: MarginalSparseSeqDistribution,
  pub msg_from_child: MarginalSparseSeqDistribution,
  pub transmission: Option<Vec<(usize, usize)>>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MarginalSparseSeqDistribution {
  /// probability vector for each variable position collecting information from children
  pub variable: BTreeMap<usize, VarPos>,

  pub variable_indel: BTreeMap<(usize, usize), Deletion>,

  /// probability vector for the state of fixed positions based on information from children
  pub fixed: BTreeMap<AsciiChar, Array1<f64>>,

  pub fixed_counts: Composition,

  /// Total log likelihood
  pub log_lh: f64,
}

impl Default for MarginalSparseSeqDistribution {
  fn default() -> Self {
    Self {
      variable: btreemap! {},
      variable_indel: btreemap! {},
      fixed: btreemap! {},
      fixed_counts: Composition::new(std::iter::empty::<AsciiChar>(), AsciiChar::from_byte_unchecked(b'-')),
      log_lh: 0.0,
    }
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FitchSeqDistribution {
  /// probability vector for each variable position collecting information from children
  pub variable: BTreeMap<usize, StateSet>,

  pub variable_indel: BTreeMap<(usize, usize), Deletion>,

  pub composition: Composition,

  /// The character chosen for each variable site during the Fitch forward pass
  pub chosen_state: BTreeMap<usize, AsciiChar>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct VarPos {
  pub dis: Array1<f64>, // array of floats of size 'alphabet'
  pub state: AsciiChar, // exact reference state for this sparse position
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
