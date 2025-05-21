use crate::representation::seq::Seq;
use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct InDel {
  pub range: (usize, usize),
  pub seq: Seq,
  pub deletion: bool, // deletion if True, insertion if False
}

impl InDel {
  pub fn del(range: (usize, usize), seq: impl Into<Seq>) -> Self {
    Self::new(range, seq, true)
  }

  pub fn ins(range: (usize, usize), seq: impl Into<Seq>) -> Self {
    Self::new(range, seq, false)
  }

  pub fn new(range: (usize, usize), seq: impl Into<Seq>, deletion: bool) -> Self {
    let seq = seq.into();
    assert!(range.0 <= range.1);
    assert_eq!(seq.len(), range.1 - range.0);
    Self { range, seq, deletion }
  }
}

impl fmt::Display for InDel {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    let delta_str = if self.deletion {
      format!("{} -> {}", &self.seq.as_str(), "-".repeat(self.seq.len()))
    } else {
      format!("{} -> {}", "-".repeat(self.seq.len()), &self.seq.as_str())
    };
    write!(f, "{}--{}: {}", self.range.0, self.range.1, delta_str)
  }
}
