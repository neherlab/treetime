use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter, EnumString};
use treetime_ops::{FftConvolve, LogScaleMultiply, NdarrayConvolve, PointwiseMultiply, RiemannConvolve};

pub use treetime_ops::traits::{ConvolveAlgo as Algo, MultiplyAlgo};

/// Available convolution algorithms for testing
#[derive(
  Debug,
  Clone,
  Copy,
  PartialEq,
  Eq,
  PartialOrd,
  Ord,
  Hash,
  Serialize,
  Deserialize,
  Display,
  EnumString,
  EnumIter,
  ValueEnum,
)]
#[serde(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
#[clap(rename_all = "kebab-case")]
pub enum ConvolutionAlgorithm {
  All,
  Riemann,
  NdarrayConv,
  NdarrayConvFft,
}

impl ConvolutionAlgorithm {
  /// Get all available algorithms (excluding the All meta-variant)
  pub fn all() -> Vec<Self> {
    Self::iter().filter(|a| *a != Self::All).collect()
  }

  /// Expand algorithm list, replacing All with actual algorithms
  pub fn expand(algorithms: &[Self]) -> Vec<Self> {
    if algorithms.contains(&Self::All) {
      Self::all()
    } else {
      algorithms.to_vec()
    }
  }

  /// Instantiate the algorithm implementation
  pub fn instantiate(&self) -> Box<dyn Algo> {
    match self {
      Self::All => panic!("Cannot instantiate All meta-variant"),
      Self::Riemann => Box::new(RiemannConvolve),
      Self::NdarrayConv => Box::new(NdarrayConvolve),
      Self::NdarrayConvFft => Box::new(FftConvolve),
    }
  }
}

/// Available multiplication algorithms for testing
#[derive(
  Debug,
  Clone,
  Copy,
  PartialEq,
  Eq,
  PartialOrd,
  Ord,
  Hash,
  Serialize,
  Deserialize,
  Display,
  EnumString,
  EnumIter,
  ValueEnum,
)]
#[serde(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
#[clap(rename_all = "kebab-case")]
pub enum MultiplicationAlgorithm {
  All,
  NaiveMultiplication,
  LogScaleMultiplication,
}

impl MultiplicationAlgorithm {
  pub fn all() -> Vec<Self> {
    Self::iter().filter(|a| *a != Self::All).collect()
  }

  pub fn expand(algorithms: &[Self]) -> Vec<Self> {
    if algorithms.contains(&Self::All) {
      Self::all()
    } else {
      algorithms.to_vec()
    }
  }

  pub fn instantiate(&self) -> Box<dyn MultiplyAlgo> {
    match self {
      Self::All => panic!("Cannot instantiate All meta-variant"),
      Self::NaiveMultiplication => Box::new(PointwiseMultiply),
      Self::LogScaleMultiplication => Box::new(LogScaleMultiply),
    }
  }
}
