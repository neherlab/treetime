use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use treetime_schema::schema_defaults;

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ZeroBasedArgs {
  /// Report sequence positions counted from 0 instead of 1.
  ///
  /// Applies to the reported mutations and positions. An input table with positions, such as `--drms`, is always
  /// counted from 1.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Input data"))]
  pub zero_based: bool,
}
