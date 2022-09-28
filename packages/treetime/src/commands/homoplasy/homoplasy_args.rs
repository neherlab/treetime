#![allow(clippy::large_enum_variant)]
#![allow(clippy::struct_excessive_bools)]

use crate::gtr::get_gtr::GtrModelName;
use clap::{ArgEnum, Parser, ValueHint};
use std::fmt::Debug;
use std::path::PathBuf;
use std::str::FromStr;

#[derive(Parser, Debug)]
pub struct TreetimeHomoplasyArgs {
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
  /// Example: treetime homoplasy seq1.fasta seq2.fasta
  #[clap(long, visible_alias("aln"))]
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(hide_long_help = true, hide_short_help = true)]
  pub aln: Option<PathBuf>,

  /// Only for vcf input: fasta file of the sequence the VCF was mapped to.
  #[clap(long, short = 'r')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub vcf_reference: Option<PathBuf>,

  /// Name of file containing the tree in newick, nexus, or phylip format.
  ///
  /// If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
  #[clap(long, short = 't')]
  #[clap(value_hint = ValueHint::FilePath)]
  pub tree: Option<PathBuf>,

  /// Number of constant sites not included in alignment
  #[clap(long = "const")]
  pub constant_sites: Option<usize>,

  /// rescale branch lengths
  #[clap(long)]
  pub rescale: bool,

  /// generate a more detailed report
  #[clap(long)]
  pub detailed: Option<String>,

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

  /// Use aminoacid alphabet
  #[clap(long)]
  pub aa: bool,

  /// Zero-based mutation indexing
  #[clap(long)]
  pub zero_based: bool,

  /// number of mutations/nodes that are printed to screen
  #[clap(long, short = 'n')]
  pub num_mut: Option<usize>,

  /// TSV file containing DRM info. Columns headers: GENOMIC_POSITION, ALT_BASE, DRUG, GENE, SUBSTITUTION
  #[clap(long)]
  #[clap(value_hint = ValueHint::FilePath)]
  pub drms: Option<PathBuf>,

  /// Directory to write the output to
  #[clap(long, short = 'O')]
  pub outdir: PathBuf,

  /// Random seed
  #[clap(long)]
  pub seed: Option<u64>,
}
