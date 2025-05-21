use crate::cli::verbosity::Verbosity;
use clap::{Parser, ValueEnum, ValueHint};
use lazy_static::lazy_static;
use std::path::{Path, PathBuf};
use treetime::io::fs::extension;
use treetime::utils::clap_styles::styles;

#[derive(Copy, Clone, Debug, ValueEnum)]
#[value(rename_all = "kebab-case")]
pub enum TreeFormat {
  Auspice,
  MatJson,
  MatPb,
  Newick,
  Nexus,
  PhyloGraph,
  Phyloxml,
  PhyloxmlJson,
}

#[derive(Parser, Debug)]
#[clap(name = "convert")]
#[clap(author, version)]
#[clap(verbatim_doc_comment)]
#[clap(styles = styles())]
/// Read and write Usher MAT files
///
/// * https://github.com/yatisht/usher
/// * https://usher-wiki.readthedocs.io/
pub struct Args {
  /// Path to input file
  ///
  /// Accepts plain or compressed files. If a compressed file is provided, it will be transparently
  /// decompressed. Supported compression formats: "gz", "bz2", "xz", "zst". Decompressor is chosen based on file
  /// extension.
  ///
  /// Omit this argument or use special value "-" to read uncompressed data from standard input (stdin).
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 1)]
  #[clap(default_value = "-")]
  pub input: PathBuf,

  /// Path to output file
  ///
  /// If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. Compressor is chosen based on file extension.
  ///
  /// Omit this argument or use special value "-" to write uncompressed data to standard output (stdout).
  #[clap(long, short = 'o')]
  #[clap(display_order = 1)]
  #[clap(value_hint = ValueHint::AnyPath)]
  #[clap(default_value = "-")]
  pub output: PathBuf,

  /// Input format to read
  ///
  /// Provide this argument if automatic format detection fails or if you want to override it.
  #[clap(long, short = 'r', value_enum)]
  #[clap(display_order = 2)]
  pub input_format: Option<TreeFormat>,

  /// Output format to write
  ///
  /// Provide this argument if automatic format detection fails or if you want to override it.
  #[clap(long, short = 'w', value_enum)]
  #[clap(display_order = 2)]
  pub output_format: Option<TreeFormat>,

  /// Make output more quiet or more verbose
  #[clap(flatten, next_help_heading = "Verbosity")]
  pub verbosity: Verbosity,
}

lazy_static! {
  static ref VERBOSITIES: &'static [&'static str] = &["off", "error", "warn", "info", "debug", "trace"];
}

pub fn guess_tree_format_from_filename(filepath: impl AsRef<Path>) -> Option<TreeFormat> {
  let filepath = filepath.as_ref();
  let ext = extension(filepath).map(|s| s.to_lowercase());
  match ext.as_deref() {
    Some("auspice.json") => Some(TreeFormat::Auspice),
    Some("graph.json") => Some(TreeFormat::PhyloGraph),
    Some("mat.json") => Some(TreeFormat::MatPb),
    Some("mat.pb") => Some(TreeFormat::MatPb),
    Some("nex | nexus") => Some(TreeFormat::Nexus),
    Some("nwk" | "newick") => Some(TreeFormat::Newick),
    Some("phylo.xml") => Some(TreeFormat::Phyloxml),
    _ => None,
  }
}
