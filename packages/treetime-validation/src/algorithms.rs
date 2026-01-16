use clap::ValueEnum;
use eyre::Report;
use serde::{Deserialize, Serialize};
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter, EnumString};
use treetime_ops::{AggressiveMultiply, FftConvolve, LogScaleMultiply, NdarrayConvolve, PointwiseMultiply, RiemannConvolve};
use treetime_utils::make_error;

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

  /// Instantiate the algorithm implementation.
  ///
  /// Returns error if called on `All` meta-variant (use `expand()` first).
  pub fn instantiate(&self) -> Result<Box<dyn Algo>, Report> {
    match self {
      Self::All => make_error!("Cannot instantiate All meta-variant; use expand() first"),
      Self::Riemann => Ok(Box::new(RiemannConvolve)),
      Self::NdarrayConv => Ok(Box::new(NdarrayConvolve)),
      Self::NdarrayConvFft => Ok(Box::new(FftConvolve)),
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
  AggressiveMultiplication,
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

  /// Instantiate the algorithm implementation.
  ///
  /// Returns error if called on `All` meta-variant (use `expand()` first).
  pub fn instantiate(&self) -> Result<Box<dyn MultiplyAlgo>, Report> {
    match self {
      Self::All => make_error!("Cannot instantiate All meta-variant; use expand() first"),
      Self::NaiveMultiplication => Ok(Box::new(PointwiseMultiply)),
      Self::LogScaleMultiplication => Ok(Box::new(LogScaleMultiply)),
      Self::AggressiveMultiplication => Ok(Box::new(AggressiveMultiply)),
    }
  }
}
