use ndarray::Array1;
use serde::{Deserialize, Serialize};
use treetime_utils::array::serde::{array1_as_vec, array1_from_vec};

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct DiscreteNodeData {
  pub observed: Option<usize>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub profile: Array1<f64>,
  pub log_lh: f64,
}

impl DiscreteNodeData {
  pub fn from_observed(index: usize, n_states: usize) -> Self {
    let mut profile = Array1::zeros(n_states);
    profile[index] = 1.0;
    Self {
      observed: Some(index),
      profile,
      log_lh: 0.0,
    }
  }

  pub fn missing(n_states: usize) -> Self {
    let profile = Array1::from_elem(n_states, 1.0 / n_states as f64);
    Self {
      observed: None,
      profile,
      log_lh: 0.0,
    }
  }
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct DiscreteEdgeData {
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub msg_to_child: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub msg_to_parent: Array1<f64>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  pub msg_from_child: Array1<f64>,
  pub log_lh: f64,
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;
  use ndarray::array;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_from_observed_creates_one_hot_profile() {
    let node = DiscreteNodeData::from_observed(1, 4);

    let expected = array![0.0, 1.0, 0.0, 0.0];
    assert_eq!(node.observed, Some(1));
    assert_abs_diff_eq!(node.profile, expected, epsilon = 1e-10);
  }

  #[test]
  fn test_missing_creates_uniform_profile() {
    let node = DiscreteNodeData::missing(4);

    let expected = array![0.25, 0.25, 0.25, 0.25];
    assert_eq!(node.observed, None);
    assert_abs_diff_eq!(node.profile, expected, epsilon = 1e-10);
  }

  #[test]
  fn test_default_node_data() {
    let node = DiscreteNodeData::default();

    assert_eq!(node.observed, None);
    assert!(node.profile.is_empty());
    assert_abs_diff_eq!(node.log_lh, 0.0, epsilon = 1e-10);
  }

  #[test]
  fn test_default_edge_data() {
    let edge = DiscreteEdgeData::default();

    assert!(edge.msg_to_child.is_empty());
    assert!(edge.msg_to_parent.is_empty());
    assert!(edge.msg_from_child.is_empty());
    assert_abs_diff_eq!(edge.log_lh, 0.0, epsilon = 1e-10);
  }
}
