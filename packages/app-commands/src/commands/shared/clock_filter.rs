use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use smart_default::SmartDefault;
use treetime_schema::schema_defaults;

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ClockFilterArgs {
  /// Ignore tips that do not follow a loose clock: the number of interquartile ranges from the root-to-tip regression
  /// beyond which a tip is an outlier.
  ///
  /// 0 switches the filter off.
  #[default = 3.0]
  #[cfg_attr(
    feature = "clap",
    clap(long, default_value_t = ClockFilterArgs::default().clock_filter, help_heading = "Molecular clock")
  )]
  pub clock_filter: f64,

  /// Allow negative clock rates.
  ///
  /// By default, rates are forced to be positive. For trees with little temporal signal, removing this restriction
  /// gives essentially mid-point rooting.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Molecular clock"))]
  pub allow_negative_rate: bool,
}
