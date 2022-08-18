#![allow(clippy::large_enum_variant)]
#![allow(clippy::struct_excessive_bools)]

use clap::{Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;
use std::str::FromStr;

#[derive(Parser, Debug)]
pub struct TreetimeMugrationArgs {
  /// Name of file containing the tree in newick, nexus, or phylip format.
  ///
  /// If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
  #[clap(long, short = 't')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub tree: Option<PathBuf>,

  /// Attribute to reconstruct, e.g. country
  #[clap(long)]
  pub attribute: Option<String>,

  /// CSV or TSV file with discrete characters. #name,country,continent taxon1,micronesia,oceania ...
  #[clap(long, short = 's')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub states: Option<PathBuf>,

  /// CSV or TSV file with probabilities of that a randomly sampled sequence at equilibrium has a particular state. E.g. population of different continents or countries. E.g.: #country,weight micronesia,0.1 ...
  #[clap(long, short = 'w')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub weights: Option<PathBuf>,

  /// Label of the column to be used as taxon name
  #[clap(long)]
  pub name_column: Option<String>,

  /// Output confidence of mugration inference
  #[clap(long)]
  #[clap(value_hint = ValueHint::AnyPath)]
  pub confidence: Option<PathBuf>,

  /// Pseudo-counts. Higher numbers results in 'flatter' models.
  #[clap(long)]
  pub pc: Option<f64>,

  /// String indicating missing data
  #[clap(long)]
  pub missing_data: Option<String>,

  /// Rough estimate of how many more events would have been observed if sequences represented an
  /// even sample. This should be roughly the (1-sum_i p_i^2)/(1-sum_i t_i^2), where p_i are the
  /// equilibrium frequencies and t_i are apparent ones.(or rather the time spent in a particular
  #[clap(long)]
  pub sampling_bias_correction: Option<String>,

  /// Directory to write the output to
  #[clap(long, short = 'O')]
  pub outdir: PathBuf,

  /// Random seed
  #[clap(long)]
  pub seed: Option<u64>,
}