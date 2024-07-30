use clap::{AppSettings, ArgEnum, Parser, ValueHint};
use clap_verbosity_flag::{Verbosity, WarnLevel};
use lazy_static::lazy_static;
use log::LevelFilter;
use std::path::{Path, PathBuf};
use treetime::io::fs::extension;

#[derive(Copy, Clone, Debug, ArgEnum)]
#[clap(rename = "kebab-case")]
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
#[clap(name = "treetime", trailing_var_arg = true)]
#[clap(author, version)]
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
#[clap(verbatim_doc_comment)]
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
  #[clap(long, short = 'r', arg_enum)]
  #[clap(display_order = 2)]
  pub input_format: Option<TreeFormat>,

  /// Output format to write
  ///
  /// Provide this argument if automatic format detection fails or if you want to override it.
  #[clap(long, short = 'w', arg_enum)]
  #[clap(display_order = 2)]
  pub output_format: Option<TreeFormat>,

  /// Make output more quiet or more verbose
  #[clap(flatten)]
  pub verbose: Verbosity<WarnLevel>,

  /// Set verbosity level
  #[clap(long, global = true, conflicts_with = "verbose", conflicts_with = "silent", possible_values(VERBOSITIES.iter()))]
  #[clap(display_order = 999)]
  pub verbosity: Option<LevelFilter>,

  /// Disable all console output. Same as --verbosity=off
  #[clap(long, global = true, conflicts_with = "verbose", conflicts_with = "verbosity")]
  #[clap(display_order = 999)]
  pub silent: bool,
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
