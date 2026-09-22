use clap::ValueEnum;
use eyre::Report;
use serde::{Deserialize, Serialize};
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter, EnumString};
use treetime_ops::{
  AggressiveMultiply, FftConvolve, LogScaleMultiply, NdarrayConvolve, PointwiseMultiply, RiemannConvolve,
};
use treetime_utils::make_error;

pub use treetime_ops::traits::{ConvolveAlgo as Algo, MultiplyAlgo};

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
  pub(crate) fn all() -> Vec<Self> {
    Self::iter().filter(|a| *a != Self::All).collect()
  }

  pub(crate) fn expand(algorithms: &[Self]) -> Vec<Self> {
    if algorithms.contains(&Self::All) {
      Self::all()
    } else {
      algorithms.to_vec()
    }
  }

  pub(crate) fn instantiate(self) -> Result<Box<dyn Algo>, Report> {
    match self {
      Self::All => make_error!("Cannot instantiate All meta-variant; use expand() first"),
      Self::Riemann => Ok(Box::new(RiemannConvolve)),
      Self::NdarrayConv => Ok(Box::new(NdarrayConvolve)),
      Self::NdarrayConvFft => Ok(Box::new(FftConvolve)),
    }
  }
}

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
  Pointwise,
  LogScale,
  Aggressive,
}

impl MultiplicationAlgorithm {
  pub(crate) fn all() -> Vec<Self> {
    Self::iter().filter(|a| *a != Self::All).collect()
  }

  pub(crate) fn expand(algorithms: &[Self]) -> Vec<Self> {
    if algorithms.contains(&Self::All) {
      Self::all()
    } else {
      algorithms.to_vec()
    }
  }

  pub(crate) fn instantiate(self) -> Result<Box<dyn MultiplyAlgo>, Report> {
    match self {
      Self::All => make_error!("Cannot instantiate All meta-variant; use expand() first"),
      Self::Pointwise => Ok(Box::new(PointwiseMultiply)),
      Self::LogScale => Ok(Box::new(LogScaleMultiply)),
      Self::Aggressive => Ok(Box::new(AggressiveMultiply)),
    }
  }
}
