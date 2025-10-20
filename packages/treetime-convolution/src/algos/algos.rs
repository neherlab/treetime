use crate::algos::ndarray_conv::ndarray_conv::NdarrayAlgo;
use crate::algos::ndarray_conv_fft::ndarray_conv_fft::NdarrayConvFftAlgo;
use crate::algos::riemann::riemann::RiemannAlgo;
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
  Riemann,
  Ndarray,
  NdarrayFft,
}

impl ConvolutionAlgorithm {
  /// Get all available algorithms
  pub fn all() -> Vec<Self> {
    Self::iter().collect()
  }

  /// Instantiate the algorithm implementation
  pub fn instantiate(&self) -> Box<dyn Algo> {
    match self {
      Self::Riemann => Box::new(RiemannAlgo),
      Self::Ndarray => Box::new(NdarrayAlgo),
      Self::NdarrayFft => Box::new(NdarrayConvFftAlgo),
    }
  }
}

pub trait Algo: Send + Sync {
  fn name(&self) -> &'static str;

  fn convolve(
    &self,
    input_grid: &Array1<f64>,
    f_values: &Array1<f64>,
    g_values: &Array1<f64>,
    output_grid: &Array1<f64>,
  ) -> Result<Array1<f64>, Report>;
}
