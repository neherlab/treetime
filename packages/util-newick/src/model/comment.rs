use crate::model::value::NewickValue;
use deser::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum NewickComment {
  Beast(Vec<(String, NewickValue)>),
  Nhx(Vec<(String, NewickValue)>),
  MrBayesMcmc(MrBayesComment),
  Plain(String),
}

impl NewickComment {
  pub fn pairs(&self) -> &[(String, NewickValue)] {
    match self {
      Self::Beast(pairs) | Self::Nhx(pairs) => pairs,
      Self::MrBayesMcmc(_) | Self::Plain(_) => &[],
    }
  }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct MrBayesComment {
  pub kind: MrBayesKind,
  pub name: String,
  pub values: Vec<String>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum MrBayesKind {
  E,
  B,
  N,
}

impl MrBayesKind {
  pub const fn letter(self) -> &'static str {
    match self {
      Self::E => "E",
      Self::B => "B",
      Self::N => "N",
    }
  }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct NodeComment {
  pub position: LabelSide,
  pub comment: NewickComment,
}

impl NodeComment {
  pub fn new(position: LabelSide, comment: NewickComment) -> Self {
    Self { position, comment }
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum LabelSide {
  BeforeLabel,
  AfterLabel,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct EdgeComment {
  pub field: EdgeField,
  pub side: ValueSide,
  pub comment: NewickComment,
}

impl EdgeComment {
  pub fn new(field: EdgeField, side: ValueSide, comment: NewickComment) -> Self {
    Self { field, side, comment }
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum EdgeField {
  Length,
  Support,
  Probability,
}

impl EdgeField {
  pub const ALL: [Self; 3] = [Self::Length, Self::Support, Self::Probability];

  pub const fn index(self) -> usize {
    match self {
      Self::Length => 0,
      Self::Support => 1,
      Self::Probability => 2,
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum ValueSide {
  BeforeValue,
  AfterValue,
}

pub(crate) fn annotation_pairs<'a>(
  comments: impl Iterator<Item = &'a NewickComment>,
) -> impl Iterator<Item = (&'a str, &'a NewickValue)> {
  comments.flat_map(|comment| comment.pairs().iter().map(|(key, value)| (key.as_str(), value)))
}
