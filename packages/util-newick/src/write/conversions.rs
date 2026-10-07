use crate::dialect::{NewickAnnotations, NewickDialect};

pub fn conversion(dialect: NewickDialect, data: DataKind) -> Conversion {
  let NewickDialect { structure, annotations } = dialect;
  let (holds, unsupported) = match data {
    DataKind::BeastComments => (annotations == NewickAnnotations::Beast, Conversion::Drop),
    DataKind::NhxComments => (annotations == NewickAnnotations::Nhx, Conversion::Drop),
    DataKind::MrBayesComments => (annotations == NewickAnnotations::MrBayes, Conversion::Drop),
    DataKind::AmpersandPlainComments => (!annotations.reserves_annotations(), Conversion::Fail),
    DataKind::RootingPlainComments => (!dialect.rooting(), Conversion::Fail),
    DataKind::WeightPlainComments => (!structure.weight(), Conversion::Fail),
    DataKind::Rooting => (dialect.rooting(), Conversion::Drop),
    DataKind::Weight => (structure.weight(), Conversion::Drop),
    DataKind::FieldSupport | DataKind::Probability | DataKind::FieldComments => {
      (structure.field_count() > 1, Conversion::Drop)
    },
    DataKind::HybridNodes => (structure.hybrid_tags(), Conversion::Fail),
  };
  if holds { Conversion::Keep } else { unsupported }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DataKind {
  BeastComments,
  NhxComments,
  MrBayesComments,
  AmpersandPlainComments,
  RootingPlainComments,
  WeightPlainComments,
  Rooting,
  Weight,
  FieldSupport,
  Probability,
  FieldComments,
  HybridNodes,
}

impl DataKind {
  pub const ALL: [Self; 12] = [
    Self::BeastComments,
    Self::NhxComments,
    Self::MrBayesComments,
    Self::AmpersandPlainComments,
    Self::RootingPlainComments,
    Self::WeightPlainComments,
    Self::Rooting,
    Self::Weight,
    Self::FieldSupport,
    Self::Probability,
    Self::FieldComments,
    Self::HybridNodes,
  ];
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Conversion {
  Keep,
  Drop,
  Fail,
}
