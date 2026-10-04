use itertools::Itertools;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_primitives::{AsciiChar, Seq};

#[derive(Clone, Debug)]
pub(crate) struct SeqOverlay {
  base: Arc<Seq>,
  overrides: BTreeMap<usize, AsciiChar>,
}

impl SeqOverlay {
  pub(crate) fn new(base: Arc<Seq>, overrides: BTreeMap<usize, AsciiChar>) -> Self {
    Self { base, overrides }
  }

  pub(crate) fn to_seq(&self) -> Seq {
    let mut seq = (*self.base).clone();
    for (&pos, &state) in &self.overrides {
      seq[pos] = state;
    }
    seq
  }

  pub(crate) fn count_differences(&self, other: &Self) -> usize {
    if Arc::ptr_eq(&self.base, &other.base) || self.base == other.base {
      return self
        .overrides
        .keys()
        .merge(other.overrides.keys())
        .dedup()
        .filter(|&&pos| self.state(pos) != other.state(pos))
        .count();
    }
    count_differing_positions(&self.to_seq(), &other.to_seq())
  }

  fn state(&self, pos: usize) -> AsciiChar {
    self.overrides.get(&pos).copied().unwrap_or_else(|| self.base[pos])
  }
}

impl From<Seq> for SeqOverlay {
  fn from(seq: Seq) -> Self {
    Self::new(Arc::new(seq), BTreeMap::new())
  }
}

fn count_differing_positions(a: &Seq, b: &Seq) -> usize {
  let shared = a.iter().zip(b.iter()).filter(|(ca, cb)| ca != cb).count();
  let length_diff = a.len().abs_diff(b.len());
  shared + length_diff
}
