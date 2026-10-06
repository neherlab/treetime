use crate::dialect::{CommentKind, NewickDialect};
use deser::{Deserialize, Serialize};

pub fn conversion(dialect: NewickDialect, data: DataKind) -> Conversion {
  let features = dialect.features();
  let (holds, unsupported) = match data {
    DataKind::BeastComments => (features.comments == CommentKind::Beast, Conversion::Drop),
    DataKind::NhxComments => (features.comments == CommentKind::Nhx, Conversion::Drop),
    DataKind::MrBayesComments => (features.comments == CommentKind::MrBayes, Conversion::Drop),
    DataKind::AmpersandPlainComments => (!features.reserves_annotations, Conversion::Fail),
    DataKind::Rooting => (features.rooting, Conversion::Drop),
    DataKind::Weight => (features.weight, Conversion::Drop),
    DataKind::FieldSupport | DataKind::Probability | DataKind::FieldComments => {
      (features.rich_fields, Conversion::Drop)
    },
    DataKind::HybridNodes => (features.hybrid_tags, Conversion::Fail),
  };
  if holds { Conversion::Keep } else { unsupported }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum DataKind {
  BeastComments,
  NhxComments,
  MrBayesComments,
  AmpersandPlainComments,
  Rooting,
  Weight,
  FieldSupport,
  Probability,
  FieldComments,
  HybridNodes,
}

impl DataKind {
  pub const ALL: [Self; 10] = [
    Self::BeastComments,
    Self::NhxComments,
    Self::MrBayesComments,
    Self::AmpersandPlainComments,
    Self::Rooting,
    Self::Weight,
    Self::FieldSupport,
    Self::Probability,
    Self::FieldComments,
    Self::HybridNodes,
  ];
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum Conversion {
  Keep,
  Drop,
  Fail,
}
