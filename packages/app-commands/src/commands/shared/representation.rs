use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct RepresentationArgs {
  /// Use the dense representation of the sequences on the tree
  ///
  /// The dense representation stores the probabilities of every state at every alignment position of each node; the
  /// sparse representation stores only the variable positions. Dense is more accurate when branches are long and many
  /// sites change, but uses more memory. When unset, the sparse representation is used.
  #[cfg_attr(feature = "clap", clap(long, help_heading = "Ancestral reconstruction"))]
  pub dense: Option<bool>,
}
