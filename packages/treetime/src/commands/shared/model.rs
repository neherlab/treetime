use crate::gtr::get_gtr::GtrModelName;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;

/// Substitution model selection shared by every command that infers or applies a rate matrix.
///
/// One flag name (`--model`, short `-g`, alias `--gtr`) replaces the earlier split between `--model`
/// (ancestral, optimize) and `--gtr` (clock, timetree). `--model` is preferred over `--gtr` because the
/// value set includes non-GTR models (for example `jtt92`). `--model-params` (alias `--gtr-params`)
/// carries model-specific `key=value` parameters.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ModelArgs {
  /// Substitution model to use
  ///
  /// `--model infer` infers a model from the data. Alternatively, specify the model type. If the
  /// specified model requires additional options, use `--model-params` to specify those.
  #[default(GtrModelName::Infer)]
  #[cfg_attr(
    feature = "clap",
    clap(long = "model", short = 'g', visible_alias = "gtr", value_enum, default_value_t = GtrModelName::Infer)
  )]
  pub model: GtrModelName,

  /// Parameters for the model selected by `--model`, given as a `key=value` list
  ///
  /// Example: `--model k80 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.
  ///
  /// See the exact definitions of the parameters in the GTR creation methods in treetime/nuc_models.py
  /// or treetime/aa_models.py
  #[cfg_attr(feature = "clap", clap(long = "model-params", visible_alias = "gtr-params"))]
  pub model_params: Vec<String>,
}
