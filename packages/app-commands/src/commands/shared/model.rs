use crate::commands::shared::mirror_enum::impl_from_mirror;
use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use smart_default::SmartDefault;
use treetime::gtr::get_gtr::GtrModelName;
use treetime::optimize::partition_rates::PartitionRates;
use treetime_schema::schema_defaults;

/// Substitution model selection shared by every command that infers or applies a rate matrix.
///
/// One flag name (`--model`, short `-g`, alias `--gtr`) serves every command; the value set includes
/// non-GTR models (for example `jtt92`). `--model-params` (alias `--gtr-params`) carries
/// model-specific `key=value` parameters.
#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ModelArgs {
  /// Substitution model to use
  ///
  /// `infer` infers a model from the data of each partition; each `--alignment` file is one partition. Alternatively,
  /// name the model type. If the named model takes parameters, use `--model-params` to set them.
  #[default(GtrModelNameCli::Infer)]
  #[cfg_attr(
    feature = "clap",
    clap(long = "model", short = 'g', visible_alias = "gtr", value_enum, default_value_t = GtrModelNameCli::Infer, help_heading = "Substitution model")
  )]
  pub model: GtrModelNameCli,

  /// Parameters for the model selected by `--model`, given as space-separated `key=value` items
  ///
  /// The parameters apply to the model of every partition; each `--alignment` file is one partition. Example: `--model
  /// k80 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.
  #[cfg_attr(
    feature = "clap",
    clap(
      long = "model-params",
      visible_alias = "gtr-params",
      num_args = 1..,
      help_heading = "Substitution model"
    )
  )]
  pub model_params: Vec<String>,
}

impl ModelArgs {
  pub fn model_name(&self) -> GtrModelName {
    self.model.into()
  }
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "GtrModelName")]
pub enum GtrModelNameCli {
  /// Infer GTR parameters from data via Fitch parsimony substitution counts.
  #[default]
  Infer,
  #[schemars(rename = "jc69")]
  #[deser(rename = "jc69")]
  JC69,
  K80,
  F81,
  #[schemars(rename = "hky85")]
  #[deser(rename = "hky85")]
  HKY85,
  T92,
  #[schemars(rename = "tn93")]
  #[deser(rename = "tn93")]
  TN93,
  #[cfg_attr(feature = "clap", value(name = "jtt92"))]
  Jtt92,
}

impl_from_mirror!(GtrModelNameCli => GtrModelName { Infer, JC69, K80, F81, HKY85, T92, TN93, Jtt92 });

/// Relative rates of the partitions, shared by the commands that optimize branch lengths.
#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct PartitionRatesArgs {
  /// Relative substitution rates of the partitions
  ///
  /// Each `--alignment` file is one partition. With two or more partitions, all partitions share the branch lengths,
  /// and each partition has one relative rate that scales them for the sites of that partition. `estimate` fits the
  /// rate of each partition with the rest of its substitution model fixed (a partition with `--model infer` starts from
  /// its fitted rate), then divides all rates by their mean weighted by the number of sites of each partition, so that
  /// this mean is 1. The branch lengths are then expected substitutions per site, averaged over all sites of all
  /// partitions. `equal` fixes the rate of every partition at 1.
  ///
  /// With one partition, this setting has no effect: the run estimates no relative rate.
  #[default(PartitionRatesCli::Estimate)]
  #[cfg_attr(
    feature = "clap",
    clap(long, value_enum, default_value_t = PartitionRatesCli::Estimate, help_heading = "Substitution model")
  )]
  pub partition_rates: PartitionRatesCli,
}

impl PartitionRatesArgs {
  pub fn partition_rates(&self) -> PartitionRates {
    self.partition_rates.into()
  }
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "PartitionRates")]
pub enum PartitionRatesCli {
  /// Estimate one relative rate per partition, normalized to a site-weighted mean of 1.
  #[default]
  Estimate,
  /// Fix the relative rate of every partition at 1.
  Equal,
}

impl_from_mirror!(PartitionRatesCli => PartitionRates { Estimate, Equal });
