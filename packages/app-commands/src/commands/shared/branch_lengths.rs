use crate::commands::shared::output_args::DivergenceUnits;
use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use treetime::optimize::indel::IndelPolicy;
use treetime_schema::schema_defaults;

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct BranchLengthArgs {
  /// Disable indel (insertion and deletion) contributions to branch-length optimization and branch-length
  /// distributions.
  ///
  /// When set, branch-length optimization uses the substitution-only likelihood, and the branch distributions of the
  /// time inference exclude the Poisson indel term. This matches standard phylogenetic tools (RAxML, IQ-TREE, PhyML,
  /// BEAST) and TreeTime 0.x. By default, indels count.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Branch lengths")
  )]
  pub no_indels: bool,

  /// Units of the divergence values in the augur node data and Auspice outputs.
  ///
  /// `mutations-per-site` (default): branch divergence as substitutions per site. `mutations`: count of reconstructed
  /// substitutions per branch that change the state: gaps, `N`, and ambiguity codes compatible with the parent state
  /// do not count. `mutations` requires reconstructed sequences, so it does not work with input branch lengths
  /// (`--branch-length-mode=input`).
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_enum, default_value_t = DivergenceUnits::default(), help_heading = "Output"))]
  pub divergence_units: DivergenceUnits,
}

pub fn indel_policy(no_indels: bool) -> IndelPolicy {
  if no_indels {
    IndelPolicy::Exclude
  } else {
    IndelPolicy::Include
  }
}
