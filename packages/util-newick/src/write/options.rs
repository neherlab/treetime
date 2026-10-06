use crate::dialect::NewickDialect;
use crate::number::NumberFormat;
use deser::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Clone, Debug, SmartDefault, Serialize, Deserialize)]
pub struct NewickWriteOptions {
  #[default(NewickDialect::Classic)]
  pub dialect: NewickDialect,
  pub branch_annotations: BranchAnnotations,
  pub support: SupportPlacement,
  #[default(true)]
  pub root_edge: bool,
  pub quoting: Quoting,
  pub spaces: Spaces,
  pub numbers: NumberFormat,
  pub indent: Option<usize>,
}

impl NewickWriteOptions {
  pub fn new(dialect: NewickDialect) -> Self {
    Self {
      dialect,
      ..Self::default()
    }
  }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum BranchAnnotations {
  #[default]
  Recorded,
  BeforeLength,
  AfterLength,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum SupportPlacement {
  #[default]
  Source,
  Label,
  Field,
  Annotation(String),
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum Quoting {
  #[default]
  WhenNeeded,
  Always,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum Spaces {
  #[default]
  Quote,
  Underscore,
}
