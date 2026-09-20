use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use treetime::gtr::get_gtr::GtrModelName;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
#[schemars(rename = "GtrModelName")]
pub enum GtrModelNameCli {
  /// Infer GTR parameters from data via Fitch parsimony substitution counts.
  #[default]
  Infer,
  #[serde(rename = "jc69")]
  JC69,
  K80,
  F81,
  #[serde(rename = "hky85")]
  HKY85,
  T92,
  #[serde(rename = "tn93")]
  TN93,
  #[cfg_attr(feature = "clap", value(name = "jtt92"))]
  Jtt92,
}

impl From<GtrModelNameCli> for GtrModelName {
  fn from(name: GtrModelNameCli) -> Self {
    match name {
      GtrModelNameCli::Infer => GtrModelName::Infer,
      GtrModelNameCli::JC69 => GtrModelName::JC69,
      GtrModelNameCli::K80 => GtrModelName::K80,
      GtrModelNameCli::F81 => GtrModelName::F81,
      GtrModelNameCli::HKY85 => GtrModelName::HKY85,
      GtrModelNameCli::T92 => GtrModelName::T92,
      GtrModelNameCli::TN93 => GtrModelName::TN93,
      GtrModelNameCli::Jtt92 => GtrModelName::Jtt92,
    }
  }
}

/// Substitution model selection shared by every command that infers or applies a rate matrix.
///
/// One flag name (`--model`, short `-g`, alias `--gtr`) replaces the earlier split between `--model`
/// (ancestral, optimize) and `--gtr` (clock, timetree). `--model` is preferred over `--gtr` because the
/// value set includes non-GTR models (for example `jtt92`). `--model-params` (alias `--gtr-params`)
/// carries model-specific `key=value` parameters.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ModelArgs {
  /// Substitution model to use
  ///
  /// `--model infer` infers a model from the data. Alternatively, specify the model type. If the
  /// specified model requires additional options, use `--model-params` to specify those.
  #[default(GtrModelNameCli::Infer)]
  #[cfg_attr(
    feature = "clap",
    clap(long = "model", short = 'g', visible_alias = "gtr", value_enum, default_value_t = GtrModelNameCli::Infer)
  )]
  pub model: GtrModelNameCli,

  /// Parameters for the model selected by `--model`, given as a `key=value` list
  ///
  /// Example: `--model k80 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.
  ///
  /// See the exact definitions of the parameters in the GTR creation methods in treetime/nuc_models.py
  /// or treetime/aa_models.py
  #[cfg_attr(feature = "clap", clap(long = "model-params", visible_alias = "gtr-params"))]
  pub model_params: Vec<String>,
}

impl ModelArgs {
  pub fn model_name(&self) -> GtrModelName {
    self.model.into()
  }
}
