use crate::algos::ndarray_conv::ndarray_conv::NdarrayAlgo;
use crate::algos::ndarray_conv_fft::ndarray_conv_fft::NdarrayConvFftAlgo;
use crate::algos::riemann::riemann::RiemannAlgo;
use crate::algos::scaled_distribution::{LogScaleMultiplicationAlgo, NaiveMultiplicationAlgo};
use clap::ValueEnum;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter, EnumString};

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
      Self::Riemann => Box::new(RiemannAlgo),
      Self::NdarrayConv => Box::new(NdarrayAlgo),
      Self::NdarrayConvFft => Box::new(NdarrayConvFftAlgo),
    }
  }
}

pub trait Algo: Send + Sync {
  fn name(&self) -> &'static str;

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report>;
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
      Self::NaiveMultiplication => Box::new(NaiveMultiplicationAlgo),
      Self::LogScaleMultiplication => Box::new(LogScaleMultiplicationAlgo),
    }
  }
}

pub trait MultiplyAlgo: Send + Sync {
  fn name(&self) -> &'static str;

  fn multiply(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report>;

  fn multiply_many(&self, dx: f64, distributions: &[&Array1<f64>]) -> Result<(Array1<f64>, f64), Report>;
}
