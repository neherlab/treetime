use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DiscreteStates {
  states: Vec<String>,
  indices: BTreeMap<String, usize>,
  missing: String,
}

impl DiscreteStates {
  #[allow(single_use_lifetimes)] // TODO: remove when anonymous lifetimes in `impl Trait` are stabilized
  pub fn from_values<'a>(values: impl Iterator<Item = &'a str>, missing: &str) -> Self {
    let mut unique: Vec<String> = values.filter(|&v| v != missing).map(|s| s.to_owned()).collect();
    unique.sort();
    unique.dedup();

    let indices = unique.iter().enumerate().map(|(i, s)| (s.clone(), i)).collect();

    Self {
      states: unique,
      indices,
      missing: missing.to_owned(),
    }
  }

  pub fn len(&self) -> usize {
    self.states.len()
  }

  pub fn is_empty(&self) -> bool {
    self.states.is_empty()
  }

  pub fn get_index(&self, name: &str) -> Option<usize> {
    if self.is_missing(name) {
      return None;
    }
    self.indices.get(name).copied()
  }

  pub fn get_name(&self, index: usize) -> &str {
    &self.states[index]
  }

  pub fn is_missing(&self, name: &str) -> bool {
    name == self.missing
  }

  pub fn iter(&self) -> impl Iterator<Item = &str> {
    self.states.iter().map(String::as_str)
  }

  pub fn missing_marker(&self) -> &str {
    &self.missing
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_from_values_sorts_and_deduplicates() {
    let values = ["USA", "Germany", "?", "China", "USA"];
    let states = DiscreteStates::from_values(values.iter().copied(), "?");

    assert_eq!(states.len(), 3);
    let names: Vec<&str> = states.iter().collect();
    assert_eq!(names, vec!["China", "Germany", "USA"]);
  }

  #[test]
  fn test_get_index_returns_none_for_missing() {
    let values = ["USA", "Germany", "?", "China"];
    let states = DiscreteStates::from_values(values.iter().copied(), "?");

    assert_eq!(states.get_index("?"), None);
  }

  #[test]
  fn test_get_index_returns_correct_index() {
    let values = ["USA", "Germany", "?", "China"];
    let states = DiscreteStates::from_values(values.iter().copied(), "?");

    assert_eq!(states.get_index("China"), Some(0));
    assert_eq!(states.get_index("Germany"), Some(1));
    assert_eq!(states.get_index("USA"), Some(2));
  }

  #[test]
  fn test_get_name_returns_correct_name() {
    let values = ["USA", "Germany", "China"];
    let states = DiscreteStates::from_values(values.iter().copied(), "?");

    assert_eq!(states.get_name(0), "China");
    assert_eq!(states.get_name(1), "Germany");
    assert_eq!(states.get_name(2), "USA");
  }

  #[test]
  fn test_is_missing() {
    let values = ["USA", "Germany"];
    let states = DiscreteStates::from_values(values.iter().copied(), "?");

    assert!(states.is_missing("?"));
    assert!(!states.is_missing("USA"));
    assert!(!states.is_missing("unknown"));
  }

  #[test]
  fn test_get_index_returns_none_for_unknown_value() {
    let values = ["USA", "Germany"];
    let states = DiscreteStates::from_values(values.iter().copied(), "?");

    assert_eq!(states.get_index("France"), None);
  }

  #[test]
  fn test_empty_values() {
    let values: [&str; 0] = [];
    let states = DiscreteStates::from_values(values.iter().copied(), "?");

    assert!(states.is_empty());
    assert_eq!(states.len(), 0);
  }

  #[test]
  fn test_missing_marker() {
    let values = ["USA"];
    let states = DiscreteStates::from_values(values.iter().copied(), "N/A");

    assert_eq!(states.missing_marker(), "N/A");
  }
}
