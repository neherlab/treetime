use crate::Seq;
use deser::{Deserialize, Serialize};

#[derive(Clone, Default, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlignmentRecord {
  pub name: String,
  pub seq: Seq,
}
