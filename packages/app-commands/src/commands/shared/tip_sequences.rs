use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use treetime_schema::schema_defaults;

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ImputeMissingDataArgs {
  /// Resolve ambiguous and unknown tip states (`N` and IUPAC codes such as `R`) to the most likely inferred state.
  ///
  /// Gaps are left as deletions, because they are inferred structure, not missing data. Changes that involve ambiguous
  /// characters on terminal branches then disappear from the branch mutations. Only defined for marginal
  /// reconstruction; a no-op with a warning under `--method-anc=parsimony`.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Ancestral reconstruction")
  )]
  pub impute_missing_data: bool,
}

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct TipSequenceArgs {
  /// Emit reconstructed leaf (tip) sequences in addition to internal nodes.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Ancestral reconstruction")
  )]
  pub include_leaves: bool,

  /// Same as `--include-leaves --impute-missing-data`: emit tip sequences with ambiguous and unknown states resolved to
  /// the most likely inferred state.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Ancestral reconstruction")
  )]
  pub reconstruct_tip_states: bool,

  /// Include branch mutations from or to the fully ambiguous state (`N` for nucleotides, `X` for amino acids).
  ///
  /// By default these mutations are omitted from the branch mutations of the Newick and Nexus annotations, the Auspice
  /// JSON, and the `muts` and `aa_muts` lists of augur node data. Other ambiguity codes, such as `K` or `R`, are always
  /// reported. The flag does not affect the MAT outputs, which store `N` as missing data and never contain these
  /// mutations.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, help_heading = "Ancestral reconstruction")
  )]
  pub report_ambiguous: bool,
}
