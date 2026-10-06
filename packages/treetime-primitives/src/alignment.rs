use crate::Seq;
use serde::{Deserialize, Serialize};

#[derive(Clone, Default, Debug, PartialEq, Eq, Serialize, Deserialize, deser::Serialize, deser::Deserialize)]
pub struct AlignmentRecord {
  pub name: String,
  pub seq: Seq,
}
