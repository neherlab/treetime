use crate::commands::shared::metadata::MetadataIdArgs;
use crate::commands::shared::output::{MugrationOutputSelection, OutputCoreArgs, TopologyOrderArgs};
#[cfg(feature = "clap")]
use clap::ValueHint;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;

#[derive(Debug, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Parser))]
pub struct TreetimeMugrationArgs {
  /// Name of file containing the tree in newick, nexus, or phylip format.
  ///
  /// If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
  #[cfg_attr(feature = "clap", clap(long, short = 't'))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub tree: Option<PathBuf>,

  /// Attribute to reconstruct, e.g. country
  #[cfg_attr(feature = "clap", clap(long))]
  pub attribute: String,

  /// CSV or TSV file with discrete characters. #name,country,continent taxon1,micronesia,oceania ...
  #[cfg_attr(feature = "clap", clap(long = "metadata", visible_alias = "states", short = 's'))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub metadata: PathBuf,

  /// CSV or TSV file with probabilities of that a randomly sampled sequence at equilibrium has a particular state. E.g. population of different continents or countries. E.g.: #country,weight micronesia,0.1 ...
  #[cfg_attr(feature = "clap", clap(long, short = 'w'))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub weights: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub metadata_id: MetadataIdArgs,

  /// Path to output state-probability-profile CSV.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, visible_alias = "confidence", value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_confidence_csv: Option<PathBuf>,

  /// Pseudo-counts. Higher numbers results in 'flatter' models. Default: 1.0.
  #[cfg_attr(feature = "clap", clap(long))]
  pub pc: Option<f64>,

  /// String indicating missing data
  #[cfg_attr(feature = "clap", clap(long, default_value = "?"))]
  #[default(_code = r#""?".to_owned()"#)]
  pub missing_data: String,

  /// Portion of attribute values that is allowed to not have weights in the weights file
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 0.5))]
  #[default = 0.5]
  pub missing_weights_threshold: f64,

  /// Number of iterations for GTR model refinement from data.
  #[cfg_attr(feature = "clap", clap(long, default_value_t = 5))]
  #[default = 5]
  pub iterations: usize,

  /// Rough estimate of how many more events would have been observed if sequences represented an
  /// even sample.
  #[cfg_attr(feature = "clap", clap(long))]
  pub sampling_bias_correction: Option<f64>,

  /// Smooth the initial equilibrium frequencies with the pseudo-count before the first
  /// reconstruction pass.
  ///
  /// Off by default (TreeTime v0 builds the initial model from raw frequencies and applies the
  /// pseudo-count only as infer_gtr regularization). Enabling this flattens the prior for the first
  /// pass; it only affects weighted models.
  #[cfg_attr(feature = "clap", clap(long))]
  pub smooth_initial_pi: bool,

  /// Exclude near-uniform root positions from the equilibrium-frequency prior.
  ///
  /// Off by default (TreeTime v0 always folds the root's most-likely state into the prior). Enabling
  /// this drops root positions whose posterior carries no phylogenetic signal, removing a
  /// state-order-dependent bias at ambiguous roots.
  #[cfg_attr(feature = "clap", clap(long))]
  pub filter_uninformative_root: bool,

  /// Path to output augur-compatible node data JSON.
  ///
  /// Contains per-node discrete trait assignments, confidence profiles, entropy,
  /// the inferred substitution model, and branch state-change labels. The output
  /// is compatible with augur export v2 --node-data for Nextstrain pipeline
  /// integration.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_augur_node_data: Option<PathBuf>,

  /// Path to output GTR model JSON.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_gtr: Option<PathBuf>,

  /// Path to output traits CSV.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_traits_csv: Option<PathBuf>,

  /// Random seed
  #[cfg_attr(feature = "clap", clap(long, visible_alias = "rng-seed"))]
  pub seed: Option<u64>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub output: OutputCoreArgs,

  /// Comma-separated list of outputs to produce with `--output-all`.
  ///
  /// Restricts which outputs `--output-all` writes. Special value `all` expands to every output
  /// available for this command. Requires `--output-all`. Per-file flags are always honored
  /// regardless of this selection.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_delimiter = ',', requires = "output_all", help_heading = "Output")
  )]
  pub output_selection: Vec<MugrationOutputSelection>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub topology_order: TopologyOrderArgs,
}
