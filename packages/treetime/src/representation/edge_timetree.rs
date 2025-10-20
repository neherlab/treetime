use crate::distribution::distribution::Distribution;
use crate::graph::edge::{GraphEdge, Weighted};
use crate::io::nwk::EdgeFromNwk;
use crate::representation::graph_ancestral::EdgeAncestral;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct EdgeTimetree {
  pub branch_length: Option<f64>,
  pub clock_length: Option<f64>,
  pub branch_length_distribution: Option<Arc<Distribution>>,
}

impl GraphEdge for EdgeTimetree {}

impl Weighted for EdgeTimetree {
  fn weight(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_weight(&mut self, weight: Option<f64>) {
    self.branch_length = weight;
  }
}

impl EdgeFromNwk for EdgeTimetree {
  fn from_nwk(branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self {
      branch_length,
      clock_length: None,
      branch_length_distribution: None,
    })
  }
}

impl From<&EdgeAncestral> for EdgeTimetree {
  fn from(edge: &EdgeAncestral) -> Self {
    Self {
      branch_length: edge.branch_length,
      clock_length: None,
      branch_length_distribution: None,
    }
  }
}

impl crate::io::nwk::EdgeToNwk for EdgeTimetree {
  fn nwk_weight(&self) -> Option<f64> {
    self.branch_length
  }
}
