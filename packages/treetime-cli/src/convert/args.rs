use crate::cli::verbosity::Verbosity;
use clap::{Parser, ValueEnum, ValueHint};
use std::path::{Path, PathBuf};
use treetime_utils::init::clap_styles::styles;

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
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
/// - https://github.com/yatisht/usher
/// - https://usher-wiki.readthedocs.io/
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

const COMPOUND_EXTENSIONS: &[(&str, TreeFormat)] = &[
  (".auspice.json", TreeFormat::Auspice),
  (".graph.json", TreeFormat::PhyloGraph),
  (".mat.json", TreeFormat::MatJson),
  (".mat.pb", TreeFormat::MatPb),
  (".phylo.xml", TreeFormat::Phyloxml),
  (".phyloxml.json", TreeFormat::PhyloxmlJson),
];

pub fn guess_tree_format_from_filename(filepath: impl AsRef<Path>) -> Option<TreeFormat> {
  let name = filepath.as_ref().file_name()?.to_str()?.to_ascii_lowercase();

  for &(suffix, format) in COMPOUND_EXTENSIONS {
    if name.ends_with(suffix) {
      return Some(format);
    }
  }

  let ext = name.rsplit('.').next()?;
  match ext {
    "nex" | "nexus" => Some(TreeFormat::Nexus),
    "nwk" | "newick" => Some(TreeFormat::Newick),
    _ => None,
  }
}
