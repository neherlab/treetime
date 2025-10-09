use serde::{Deserialize, Serialize};
use std::fmt;

/// Available convolution algorithms for testing
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConvolutionAlgorithm {
  #[serde(rename = "riemann")]
  Riemann,
  #[serde(rename = "ndarray")]
  NdArray,
}

impl fmt::Display for ConvolutionAlgorithm {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    match self {
      ConvolutionAlgorithm::Riemann => write!(f, "riemann"),
      ConvolutionAlgorithm::NdArray => write!(f, "ndarray"),
    }
  }
}

impl ConvolutionAlgorithm {
  /// Get all available algorithms
  pub fn all() -> Vec<Self> {
    vec![Self::Riemann, Self::NdArray]
  }

  /// Parse from string
  pub fn from_str(s: &str) -> Result<Self, String> {
    match s.trim() {
      "riemann" => Ok(Self::Riemann),
      "ndarray" => Ok(Self::NdArray),
      other => Err(format!("Unknown algorithm: {}", other)),
    }
  }

  /// Parse multiple algorithms from comma-separated string
  pub fn parse_list(s: &str) -> Result<Vec<Self>, eyre::Report> {
    s.split(',')
      .map(|s| Self::from_str(s.trim()))
      .collect::<Result<Vec<_>, _>>()
      .map_err(|e| eyre::eyre!("{}", e))
  }
}
