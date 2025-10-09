use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter, EnumString};

/// Available convolution algorithms for testing
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Display, EnumString, EnumIter, ValueEnum)]
#[serde(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
#[clap(rename_all = "kebab-case")]
pub enum ConvolutionAlgorithm {
  Riemann,
  Ndarray,
}

impl ConvolutionAlgorithm {
  /// Get all available algorithms
  pub fn all() -> Vec<Self> {
    Self::iter().collect()
  }
}
