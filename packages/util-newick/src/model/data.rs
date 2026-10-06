use crate::model::comment::{EdgeComment, NodeComment, annotation_pairs};
use crate::model::value::NewickValue;
use deser::{Deserialize, Serialize};

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct NewickNodeData {
  name: Option<String>,
  hybrid: Option<NewickHybrid>,
  comments: Vec<NodeComment>,
}

impl NewickNodeData {
  pub fn new() -> Self {
    Self::default()
  }

  #[must_use]
  pub fn with_name(mut self, name: impl Into<String>) -> Self {
    self.name = Some(name.into());
    self
  }

  #[must_use]
  pub fn with_hybrid(mut self, hybrid: NewickHybrid) -> Self {
    self.hybrid = Some(hybrid);
    self
  }

  #[must_use]
  pub fn with_comment(mut self, comment: NodeComment) -> Self {
    self.comments.push(comment);
    self
  }

  pub fn name(&self) -> Option<&str> {
    self.name.as_deref()
  }

  pub fn set_name(&mut self, name: Option<String>) {
    self.name = name;
  }

  pub fn hybrid(&self) -> Option<&NewickHybrid> {
    self.hybrid.as_ref()
  }

  pub fn comments(&self) -> &[NodeComment] {
    &self.comments
  }

  pub fn comments_mut(&mut self) -> &mut Vec<NodeComment> {
    &mut self.comments
  }

  pub fn annotations(&self) -> impl Iterator<Item = (&str, &NewickValue)> {
    annotation_pairs(self.comments.iter().map(|comment| &comment.comment))
  }

  pub fn annotation(&self, key: &str) -> Option<&NewickValue> {
    last_value(self.annotations(), key)
  }

  pub(crate) fn into_comments(self) -> Vec<NodeComment> {
    self.comments
  }

  pub(crate) fn has_label(&self) -> bool {
    self.name.is_some() || self.hybrid.is_some()
  }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct NewickHybrid {
  pub kind: Option<String>,
  pub index: u32,
}

impl NewickHybrid {
  pub fn new(kind: Option<String>, index: u32) -> Self {
    Self { kind, index }
  }

  pub(crate) fn tag(&self, is_acceptor: bool) -> String {
    let marker = if is_acceptor { "##" } else { "#" };
    format!("{marker}{}{}", self.kind.as_deref().unwrap_or(""), self.index)
  }
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct NewickEdgeData {
  branch_length: Option<f64>,
  support: Vec<f64>,
  support_source: SupportSource,
  probability: Option<f64>,
  is_acceptor: bool,
  comments: Vec<EdgeComment>,
}

impl NewickEdgeData {
  pub fn new() -> Self {
    Self::default()
  }

  #[must_use]
  pub fn with_length(mut self, length: f64) -> Self {
    self.branch_length = Some(length);
    self
  }

  #[must_use]
  pub fn with_support(mut self, support: Vec<f64>, source: SupportSource) -> Self {
    self.support = support;
    self.support_source = source;
    self
  }

  #[must_use]
  pub fn with_probability(mut self, probability: f64) -> Self {
    self.probability = Some(probability);
    self
  }

  #[must_use]
  pub fn with_acceptor(mut self, is_acceptor: bool) -> Self {
    self.is_acceptor = is_acceptor;
    self
  }

  #[must_use]
  pub fn with_comment(mut self, comment: EdgeComment) -> Self {
    self.comments.push(comment);
    self
  }

  pub fn branch_length(&self) -> Option<f64> {
    self.branch_length
  }

  pub fn set_branch_length(&mut self, length: Option<f64>) {
    self.branch_length = length;
  }

  pub fn support(&self) -> &[f64] {
    &self.support
  }

  pub fn support_source(&self) -> SupportSource {
    self.support_source
  }

  pub fn set_support(&mut self, support: Vec<f64>, source: SupportSource) {
    self.support = support;
    self.support_source = source;
  }

  pub fn probability(&self) -> Option<f64> {
    self.probability
  }

  pub fn set_probability(&mut self, probability: Option<f64>) {
    self.probability = probability;
  }

  pub fn is_acceptor(&self) -> bool {
    self.is_acceptor
  }

  pub fn comments(&self) -> &[EdgeComment] {
    &self.comments
  }

  pub fn comments_mut(&mut self) -> &mut Vec<EdgeComment> {
    &mut self.comments
  }

  pub fn annotations(&self) -> impl Iterator<Item = (&str, &NewickValue)> {
    annotation_pairs(self.comments.iter().map(|comment| &comment.comment))
  }

  pub fn annotation(&self, key: &str) -> Option<&NewickValue> {
    last_value(self.annotations(), key)
  }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum SupportSource {
  #[default]
  Label,
  Field,
}

fn last_value<'a>(pairs: impl Iterator<Item = (&'a str, &'a NewickValue)>, key: &str) -> Option<&'a NewickValue> {
  pairs.filter(|(name, _)| *name == key).last().map(|(_, value)| value)
}
