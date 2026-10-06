use crate::dialect::NewickDialect;
use crate::model::graph::NewickGraph;
use crate::read::error::NewickWarning;
use deser::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Clone, Debug, SmartDefault, Serialize, Deserialize)]
pub struct NewickReadOptions {
  #[default(vec![NewickDialect::Classic])]
  pub dialects: Vec<NewickDialect>,
  pub mode: ReadMode,
  pub internal_label: InternalLabel,
  pub underscores_as_spaces: bool,
}

impl NewickReadOptions {
  pub fn all_dialects() -> Self {
    Self {
      dialects: NewickDialect::ALL.to_vec(),
      ..Self::default()
    }
  }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum ReadMode {
  #[default]
  Strict,
  Tolerant,
}

impl ReadMode {
  pub const fn name(self) -> &'static str {
    match self {
      Self::Strict => "strict",
      Self::Tolerant => "tolerant",
    }
  }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum InternalLabel {
  #[default]
  Auto,
  Name,
  Support,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct NewickTree {
  pub graph: NewickGraph,
  pub dialect: NewickDialect,
  pub warnings: Vec<NewickWarning>,
}
