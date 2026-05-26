use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use std::path::PathBuf;
#[cfg(feature = "clap")]
use clap::ValueHint;

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
  #[cfg_attr(feature = "clap", clap(long, short = 's'))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub states: PathBuf,

  /// CSV or TSV file with probabilities of that a randomly sampled sequence at equilibrium has a particular state. E.g. population of different continents or countries. E.g.: #country,weight micronesia,0.1 ...
  #[cfg_attr(feature = "clap", clap(long, short = 'w'))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::FilePath))]
  pub weights: Option<PathBuf>,

  /// Label of the column to be used as taxon name
  #[cfg_attr(feature = "clap", clap(long))]
  pub name_column: Option<String>,

  /// Output confidence of mugration inference
  #[cfg_attr(feature = "clap", clap(long))]
  #[cfg_attr(feature = "clap", clap(value_hint = ValueHint::AnyPath))]
  pub confidence: Option<PathBuf>,

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

  /// Directory to write the output to
  #[cfg_attr(feature = "clap", clap(long, short = 'O'))]
  pub outdir: PathBuf,
}
