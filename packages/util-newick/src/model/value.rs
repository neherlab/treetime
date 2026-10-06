use deser::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::hash::{Hash, Hasher};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[deser(rename_all = "kebab-case")]
pub enum NewickValue {
  Boolean(bool),
  Number(f64),
  NumberText(String),
  String(String),
  Color([u8; 3]),
  Array(NewickArray),
}

impl NewickValue {
  pub fn as_f64(&self) -> Option<f64> {
    match self {
      Self::Number(number) => Some(*number),
      Self::NumberText(text) => text.parse().ok(),
      Self::Boolean(_) | Self::String(_) | Self::Color(_) | Self::Array(_) => None,
    }
  }

  pub fn as_str(&self) -> Option<&str> {
    match self {
      Self::String(text) => Some(text),
      Self::Boolean(_) | Self::Number(_) | Self::NumberText(_) | Self::Color(_) | Self::Array(_) => None,
    }
  }

  fn rank(&self) -> u8 {
    match self {
      Self::Boolean(_) => 0,
      Self::Number(_) => 1,
      Self::NumberText(_) => 2,
      Self::String(_) => 3,
      Self::Color(_) => 4,
      Self::Array(_) => 5,
    }
  }
}

impl PartialEq for NewickValue {
  fn eq(&self, other: &Self) -> bool {
    self.cmp(other) == Ordering::Equal
  }
}

impl Eq for NewickValue {}

impl PartialOrd for NewickValue {
  fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
    Some(self.cmp(other))
  }
}

impl Ord for NewickValue {
  fn cmp(&self, other: &Self) -> Ordering {
    match (self, other) {
      (Self::Boolean(a), Self::Boolean(b)) => a.cmp(b),
      (Self::Number(a), Self::Number(b)) => a.to_bits().cmp(&b.to_bits()),
      (Self::NumberText(a), Self::NumberText(b)) | (Self::String(a), Self::String(b)) => a.cmp(b),
      (Self::Color(a), Self::Color(b)) => a.cmp(b),
      (Self::Array(a), Self::Array(b)) => a.cmp(b),
      _ => self.rank().cmp(&other.rank()),
    }
  }
}

impl Hash for NewickValue {
  fn hash<H: Hasher>(&self, state: &mut H) {
    self.rank().hash(state);
    match self {
      Self::Boolean(value) => value.hash(state),
      Self::Number(value) => value.to_bits().hash(state),
      Self::NumberText(value) | Self::String(value) => value.hash(state),
      Self::Color(value) => value.hash(state),
      Self::Array(values) => values.hash(state),
    }
  }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[deser(transparent)]
pub struct NewickArray(Vec<NewickValue>);

impl NewickArray {
  pub fn new(values: Vec<NewickValue>) -> Self {
    Self(values)
  }

  pub fn as_slice(&self) -> &[NewickValue] {
    &self.0
  }

  pub fn push(&mut self, value: NewickValue) {
    self.0.push(value);
  }

  pub fn into_vec(mut self) -> Vec<NewickValue> {
    std::mem::take(&mut self.0)
  }
}

impl From<Vec<NewickValue>> for NewickArray {
  fn from(values: Vec<NewickValue>) -> Self {
    Self(values)
  }
}

impl Drop for NewickArray {
  fn drop(&mut self) {
    let mut pending = std::mem::take(&mut self.0);
    while let Some(value) = pending.pop() {
      if let NewickValue::Array(mut inner) = value {
        pending.append(&mut inner.0);
      }
    }
  }
}
