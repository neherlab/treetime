use crate::substitution_model::NUCLEOTIDES;
use eyre::Report;
use itertools::Itertools;
use rand::{Rng, RngExt};
use std::ops::RangeInclusive;
use std::sync::LazyLock;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::NodeMap;
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::{make_error, make_report};

pub(crate) const UNKNOWN: AsciiChar = AsciiChar::from_byte_unchecked(b'N');

pub fn observe_tips(
  graph: &Graph,
  sequences: &NodeMap<Seq>,
  observation: &TipObservation,
  rng: &mut impl Rng,
) -> Result<NodeMap<Seq>, Report> {
  observation.validate()?;
  graph
    .get_leaves()
    .map(|leaf| {
      let key = leaf.key();
      let seq = sequences
        .get(key)
        .ok_or_else(|| make_report!("Leaf {key} has no sequence"))?;
      Ok((key, observation.observe(seq, rng)?))
    })
    .collect()
}

#[derive(Clone, Debug, PartialEq)]
pub struct TipObservation {
  pub ambiguity_rate: f64,
  pub unknown_run_rate: f64,
  pub unknown_run_lengths: RangeInclusive<usize>,
}

impl TipObservation {
  fn validate(&self) -> Result<(), Report> {
    for (what, rate) in [
      ("Ambiguity rate", self.ambiguity_rate),
      ("Rate of runs of N", self.unknown_run_rate),
    ] {
      if !(0.0..=1.0).contains(&rate) {
        return make_error!("{what} must be in [0, 1], got {rate}");
      }
    }
    if *self.unknown_run_lengths.start() == 0 || self.unknown_run_lengths.is_empty() {
      return make_error!(
        "Lengths of runs of N must be a non-empty range of positive lengths, got {:?}",
        self.unknown_run_lengths
      );
    }
    Ok(())
  }

  fn observe(&self, seq: &Seq, rng: &mut impl Rng) -> Result<Seq, Report> {
    let mut observed = seq
      .iter()
      .map(|&state| match codes_containing(state) {
        Some(codes) if rng.random_bool(self.ambiguity_rate) => codes[rng.random_range(0..codes.len())],
        _ => state,
      })
      .collect::<Seq>();
    if rng.random_bool(self.unknown_run_rate) {
      let max_length = *self.unknown_run_lengths.end();
      if max_length > observed.len() {
        return make_error!(
          "Runs of N up to {max_length} sites do not fit into a sequence of {} sites",
          observed.len()
        );
      }
      let length = rng.random_range(self.unknown_run_lengths.clone());
      let start = rng.random_range(0..=observed.len() - length);
      observed.as_mut_slice()[start..start + length].fill(UNKNOWN);
    }
    Ok(observed)
  }
}

pub(crate) fn codes_containing(state: AsciiChar) -> Option<&'static [AsciiChar]> {
  NUCLEOTIDES
    .iter()
    .position(|&nucleotide| nucleotide == state)
    .map(|index| CODES_CONTAINING[index].as_slice())
}

const AMBIGUITY_CODES: [(u8, &[u8]); 10] = [
  (b'R', b"AG"),
  (b'Y', b"CT"),
  (b'S', b"CG"),
  (b'W', b"AT"),
  (b'K', b"GT"),
  (b'M', b"AC"),
  (b'B', b"CGT"),
  (b'D', b"AGT"),
  (b'H', b"ACT"),
  (b'V', b"ACG"),
];

static CODES_CONTAINING: LazyLock<[Vec<AsciiChar>; 4]> = LazyLock::new(|| {
  NUCLEOTIDES.map(|nucleotide| {
    AMBIGUITY_CODES
      .iter()
      .filter(|(_, states)| states.contains(&nucleotide.inner()))
      .map(|&(code, _)| AsciiChar::from_byte_unchecked(code))
      .collect_vec()
  })
});
