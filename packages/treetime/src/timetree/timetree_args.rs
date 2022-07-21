#![allow(clippy::large_enum_variant)]
#![allow(clippy::struct_excessive_bools)]

use crate::ancestral::anc_args::MethodAncestral;
use crate::gtr::get_gtr::GtrModelName;
use clap::{ArgEnum, Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;
use std::str::FromStr;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum)]
#[clap(rename = "kebab-case")]
pub enum BranchLengthMode {
  Auto,
  Input,
  Joint,
  Marginal,
}

impl Default for BranchLengthMode {
  fn default() -> Self {
    Self::Auto
  }
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum)]
#[clap(rename = "kebab-case")]
pub enum TimeMarginalMode {
  Never,
  Always,
  OnlyFinal,
}

impl Default for TimeMarginalMode {
  fn default() -> Self {
    Self::Never
  }
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum)]
#[clap(rename = "kebab-case")]
pub enum RerootMode {
  LeastSquares,
  MinDev,
  Oldest,
}

impl Default for RerootMode {
  fn default() -> Self {
    Self::LeastSquares
  }
}

#[derive(Parser, Debug)]
pub struct TreetimeTimetreeArgs {
  /// Path to one or multiple FASTA files with aligned input sequences
  ///
  /// Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension. If there's multiple input files, then different files can have different compression formats.
  ///
  /// If no input files provided, the plain fasta input is read from standard input (stdin).
  ///
  /// See: https://en.wikipedia.org/wiki/FASTA_format
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 1)]
  pub input_fastas: Vec<PathBuf>,

  /// REMOVED. Use positional arguments instead.
  ///
  /// Example: treetime timetree seq1.fasta seq2.fasta
  #[clap(long, visible_alias("aln"))]
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(hide_long_help = true, hide_short_help = true)]
  pub aln: Option<PathBuf>,

  /// Name of file containing the tree in newick, nexus, or phylip format.
  ///
  /// If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
  #[clap(long, short = 't')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub tree: Option<PathBuf>,

  /// Only for vcf input: fasta file of the sequence the VCF was mapped to.
  #[clap(long, short = 'r')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub vcf_reference: Option<PathBuf>,

  /// CSV file with dates for nodes with 'node_name, date' where date is float (as in 2012.15)
  #[clap(long, short = 'd')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub dates: Option<PathBuf>,

  /// Label of the column to be used as taxon name
  #[clap(long)]
  pub name_column: Option<String>,

  /// Label of the column to be used as sampling date
  #[clap(long)]
  pub date_column: Option<String>,

  /// Length of the sequence, used to calculate expected variation in branch length. Not required if alignment is provided.
  #[clap(long)]
  pub sequence_length: Option<usize>,

  /// If specified, the rate of the molecular clock won't be optimized.
  #[clap(long)]
  pub clock_rate: Option<f64>,

  /// Standard deviation of the provided clock rate estimate
  #[clap(long)]
  pub clock_std_dev: Option<f64>,

  /// If set to 'input', the provided branch length will be used without modification. Note that branch lengths optimized by treetime are only accurate at short evolutionary distances.
  #[clap(long, arg_enum, default_value_t = BranchLengthMode::default())]
  pub branch_length_mode: BranchLengthMode,

  /// For 'false' or 'never', TreeTime uses the jointly most likely values for the divergence times.
  ///  For 'true' and 'always', it uses the marginal inference mode at every round of optimization,
  ///  for 'only-final' (or 'assign' for compatibility with previous versions) only uses the marginal
  ///  distribution in the final round.
  #[clap(long, arg_enum, default_value_t = TimeMarginalMode::default())]
  pub time_marginal: TimeMarginalMode,

  /// estimate confidence intervals of divergence times using the marginal posterior distribution,
  /// if `--time-marginal` is False (default) inferred divergence times will still be calculated
  /// using the jointly most likely tree configuration.
  #[clap(long)]
  pub confidence: bool,

  /// Don't resolve polytomies using temporal information.
  #[clap(long)]
  pub keep_polytomies: bool,

  /// use an autocorrelated molecular clock. Strength of the gaussian priors on branch specific rate
  /// deviation and the coupling of parent and offspring rates can be specified e.g. as --relax 1.0
  /// 0.5. Values around 1.0 correspond to weak priors, larger values constrain rate deviations more
  /// strongly. Coupling 0 (--relax 1.0 0) corresponds to an un-correlated clock.
  #[clap(long)]
  pub relax: Vec<f64>,

  /// maximal number of iterations the inference cycle is run. Note that for polytomy resolution and
  /// coalescence models max_iter should be at least 2
  #[clap(long)]
  pub max_iter: Option<usize>,

  /// coalescent time scale -- sensible values are on the order of the average hamming distance of
  /// contemporaneous sequences. In addition, 'opt' 'skyline' are valid options and estimate a
  /// constant coalescent rate or a piecewise linear coalescent rate history
  #[clap(long)]
  pub coalescent: Option<String>,

  /// number of grid points in skyline coalescent model
  #[clap(long)]
  pub n_skyline: Option<usize>,

  /// add posterior LH to coalescent model: use the posterior probability distributions of
  /// divergence times for estimating the number of branches when calculating the coalescent
  /// mergerrate or use inferred time before present (default).
  #[clap(long)]
  pub n_branches_posterior: Option<usize>,

  /// filename to save the plot to. Suffix will determine format (choices pdf, png, svg,
  /// default=pdf)
  #[clap(long)]
  #[clap(value_hint = ValueHint::FilePath)]
  pub plot_tree: Option<usize>,

  /// filename to save the plot to. Suffix will determine format (choices pdf, png, svg,
  /// default=pdf)
  #[clap(long)]
  #[clap(value_hint = ValueHint::FilePath)]
  pub plot_rtt: Option<usize>,

  /// add tip labels (default for small trees with <30 leaves)
  #[clap(long)]
  pub tip_labels: bool,

  /// don't show tip labels (default for trees with >=30 leaves)
  #[clap(long)]
  pub no_tip_labels: bool,

  /// ignore tips that don't follow a loose clock, 'clock-filter=number of inter-quartile ranges from
  /// regression'. Default=3.0, set to 0 to switch off.
  #[clap(long, default_value = "3.0")]
  pub clock_filter: f64,

  /// Reroot the tree using root-to-tip regression. Valid choices are 'min_dev', 'least-squares',
  /// and 'oldest'. 'least-squares' adjusts the root to minimize residuals of the root-to-tip vs
  /// sampling time regression, 'min_dev' minimizes variance of root-to-tip distances. 'least-
  /// squares' can be combined with --covariation to account for shared ancestry. Alternatively, you
  /// can specify a node name or a list of node names to be used as outgroup or use 'oldest' to
  /// reroot to the oldest node. By default, TreeTime will reroot using 'least-squares'. Use --keep-
  /// root to keep the current root.
  #[clap(long, arg_enum, default_value_t = RerootMode::default())]
  pub reroot: RerootMode,

  /// don't reroot the tree. Otherwise, reroot to minimize the the residual of the regression of
  /// root-to-tip distance and sampling time
  #[clap(long)]
  pub keep_root: bool,

  /// excess variance associated with terminal nodes accounting for overdispersion of the molecular
  /// clock
  #[clap(long)]
  pub tip_slack: Option<f64>,

  /// Account for covariation when estimating rates or rerooting using root-to-tip regression
  #[clap(long)]
  pub covariation: bool,

  /// GTR model to use
  ///
  /// '--gtr infer' will infer a model from the data. Alternatively, specify the model type. If the specified model requires additional options, use '--gtr-params' to specify those.
  #[clap(long, short = 'g', arg_enum, default_value_t = GtrModelName::default())]
  pub gtr: GtrModelName,

  /// GTR parameters for the model specified by the --gtr argument. The parameters should be feed as 'key=value' list of parameters.
  ///
  /// Example: '--gtr K80 --gtr-params kappa=0.2 pis=0.25,0.25,0.25,0.25'.
  ///
  /// See the exact definitions of the parameters in the GTR creation methods in treetime/nuc_models.py or treetime/aa_models.py
  #[clap(long)]
  pub gtr_params: Vec<String>,

  /// Method used for reconstructing ancestral sequences
  #[clap(long, arg_enum, default_value_t = MethodAncestral::default())]
  pub method_anc: MethodAncestral,

  /// Use aminoacid alphabet
  #[clap(long)]
  pub aa: bool,

  /// Do not fill terminal gaps
  #[clap(long)]
  pub keep_overhangs: bool,

  /// Zero-based mutation indexing
  #[clap(long)]
  pub zero_based: bool,

  /// Overwrite ambiguous states on tips with the most likely inferred state
  #[clap(long)]
  pub reconstruct_tip_states: bool,

  /// Include transitions involving ambiguous states
  #[clap(long)]
  pub report_ambiguous: bool,

  /// Directory to write the output to
  #[clap(long, short = 'O')]
  pub outdir: PathBuf,

  /// Random seed
  #[clap(long)]
  pub seed: Option<u64>,
}
