use crate::commands::ancestral::anc_args::MethodAncestral;
use crate::commands::timetree::timetree_args::{BranchLengthMode, RerootMode};
use crate::gtr::get_gtr::GtrModelName;
use clap::{Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;

#[derive(Parser, Debug)]
pub struct TreetimeClockArgs {
  /// Path to one or multiple FASTA files with aligned input sequences
  ///
  /// Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently
  /// decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file
  /// extension. If there's multiple input files, then different files can have different compression formats.
  ///
  /// If no input files provided, the plain fasta input is read from standard input (stdin).
  ///
  /// See: https://en.wikipedia.org/wiki/FASTA_format
  #[clap(long)]
  #[clap(value_hint = ValueHint::FilePath)]
  pub aln: Vec<PathBuf>,

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
  pub dates: PathBuf,

  /// Label of the column to be used as taxon name
  #[clap(long)]
  pub name_column: Option<String>,

  /// Label of the column to be used as sampling date
  #[clap(long)]
  pub date_column: Option<String>,

  /// Length of the sequence, used to calculate expected variation in branch length. Not required if alignment is provided.
  #[clap(long)]
  pub sequence_length: Option<usize>,

  /// GTR model to use
  ///
  /// '--gtr infer' will infer a model from the data. Alternatively, specify the model type. If the specified model requires additional options, use '--gtr-params' to specify those.
  #[clap(long, short = 'g', value_enum, default_value_t = GtrModelName::default())]
  pub gtr: GtrModelName,

  /// GTR parameters for the model specified by the --gtr argument. The parameters should be feed as 'key=value' list of parameters.
  ///
  /// Example: '--gtr K80 --gtr-params kappa=0.2 pis=0.25,0.25,0.25,0.25'.
  ///
  /// See the exact definitions of the parameters in the GTR creation methods in treetime/nuc_models.py or treetime/aa_models.py
  #[clap(long)]
  pub gtr_params: Vec<String>,

  /// If set to 'input', the provided branch length will be used without modification. Note that branch lengths optimized by treetime are only accurate at short evolutionary distances.
  #[clap(long, value_enum, default_value_t = BranchLengthMode::default())]
  pub branch_length_mode: BranchLengthMode,

  /// Method used for reconstructing ancestral sequences
  #[clap(long, value_enum, default_value_t = MethodAncestral::default())]
  pub method_anc: MethodAncestral,

  /// ignore tips that don't follow a loose clock, 'clock-filter=number of interquartile ranges from regression'. Default=3.0, set to 0 to switch off.
  #[clap(long, default_value = "3.0")]
  pub clock_filter: f64,

  /// Reroot the tree using root-to-tip regression. Valid choices are 'min_dev', 'least-squares',
  /// and 'oldest'. 'least-squares' adjusts the root to minimize residuals of the root-to-tip vs
  /// sampling time regression, 'min_dev' minimizes variance of root-to-tip distances. 'least-
  /// squares' can be combined with --covariation to account for shared ancestry. Alternatively, you
  /// can specify a node name or a list of node names to be used as outgroup or use 'oldest' to
  /// reroot to the oldest node. By default, TreeTime will reroot using 'least-squares'. Use --keep-
  /// root to keep the current root.
  #[clap(long, value_enum, default_value_t = RerootMode::default())]
  pub reroot: RerootMode,

  /// don't reroot the tree. Otherwise, reroot to minimize the the residual of the regression of
  /// root-to-tip distance and sampling time
  #[clap(long)]
  pub keep_root: bool,

  #[clap(long)]
  pub prune_short: bool,

  /// excess variance associated with terminal nodes accounting for overdispersion of the molecular
  /// clock
  #[clap(long)]
  pub tip_slack: Option<f64>,

  /// Account for covariation when estimating rates or rerooting using root-to-tip regression
  #[clap(long)]
  pub covariation: bool,

  /// By default, rates are forced to be positive. For trees with little temporal signal it is advisable to remove this restriction to achieve essentially mid-point rooting.
  #[clap(long)]
  pub allow_negative_rate: bool,

  /// filename to save the plot to. Suffix will determine format (choices pdf, png, svg,
  /// default=pdf)
  #[clap(long)]
  #[clap(value_hint = ValueHint::FilePath)]
  pub plot_rtt: Option<usize>,

  /// Directory to write the output to
  #[clap(long, short = 'O')]
  pub outdir: PathBuf,

  /// Random seed
  #[clap(long)]
  pub seed: Option<u64>,
}
