use crate::seq::serde::{serde_deserialize_seq, serde_serialize_seq};
use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct InDel {
  pub range: (usize, usize),
  #[serde(serialize_with = "serde_serialize_seq", deserialize_with = "serde_deserialize_seq")]
  pub seq: Vec<char>,
  pub deletion: bool, // deletion if True, insertion if False
}

impl fmt::Display for InDel {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    let delta_str = if self.deletion {
      format!("{} -> {}", String::from_iter(&self.seq), "-".repeat(self.seq.len()))
    } else {
      format!("{} -> {}", "-".repeat(self.seq.len()), String::from_iter(&self.seq))
    };
    write!(f, "{}--{}: {}", self.range.0, self.range.1, delta_str)
  }
}
