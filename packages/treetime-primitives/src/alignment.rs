use crate::Seq;

#[derive(Clone, Default, Debug, PartialEq, Eq)]
pub struct AlignmentRecord {
  pub name: String,
  pub seq: Seq,
}
