use crate::distribution::reference::convolution_test::functions::exponential::ExponentialConvInput;
use crate::distribution::reference::convolution_test::functions::gaussian::GaussianConvInput;
use crate::distribution::reference::convolution_test::traits::ConvInput;
use clap::ValueEnum;
use eyre::Report;
use serde::{Deserialize, Serialize};
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter};

#[derive(Copy, Clone, Debug, Default, Display, ValueEnum, Serialize, Deserialize, EnumIter)]
#[serde(rename_all = "kebab-case")]
#[clap(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum FunctionType {
  #[default]
  Gaussian,
  Exponential,
}

impl FunctionType {
  /// Get all available function types
  pub fn all() -> Vec<Self> {
    Self::iter().collect()
  }
}

// Callers should match on FunctionType directly:
// match function_type {
//   FunctionType::Gaussian => GaussianConvInput::new(test_cases)?,
//   FunctionType::Exponential => ExponentialConvInput::new(test_cases)?,
// }
