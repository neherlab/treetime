use crate::Seq;
use serde::{Deserialize, Serialize};

/// One alignment sequence with its sample name: the slim, domain-facing form the reconstruction
/// pipeline consumes, holding only the name the tree matches on and the sequence itself.
#[derive(Clone, Default, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlignmentRecord {
  pub name: String,
  pub seq: Seq,
}
