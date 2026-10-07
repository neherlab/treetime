use crate::dialect::NewickDialect;
use crate::model::graph::NewickGraph;
use crate::read::error::NewickWarning;
use smart_default::SmartDefault;

#[derive(Clone, Debug, SmartDefault)]
pub struct NewickReadOptions {
  pub dialect: NewickDialect,
  pub mode: ReadMode,
  pub internal_label: InternalLabel,
  pub underscores_as_spaces: bool,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
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

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum InternalLabel {
  #[default]
  Auto,
  Name,
  Support,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NewickTree {
  pub graph: NewickGraph,
  pub warnings: Vec<NewickWarning>,
}
