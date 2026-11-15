use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::input_warnings::warn_duplicate_names;
#[cfg(feature = "clap")]
use clap::ValueHint;
#[cfg(feature = "clap")]
use clap::builder::PossibleValue;
use derive_more::{Display, FromStr};
use deser::adapters::DisplayFromStr;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::{JsonSchema, Schema, SchemaGenerator, json_schema};
use smart_default::SmartDefault;
use std::borrow::Cow;
use std::path::{Path, PathBuf};
use std::sync::LazyLock;
use treetime::progress::{LogSink, RunWarningKind};
use treetime_io::nwk::NewickDialect;
use treetime_io::nwk::{NwkParse, TREE_DIALECT_DEFAULT};
use treetime_io::tree::tree_read_file;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

static TREE_DIALECTS: LazyLock<Vec<TreeDialectArg>> =
  LazyLock::new(|| NewickDialect::pairs().map(TreeDialectArg).collect());

static TREE_DIALECT_NAMES: LazyLock<Vec<String>> =
  LazyLock::new(|| TREE_DIALECTS.iter().map(ToString::to_string).collect());

pub fn read_input_tree(path: &Path, dialect: NewickDialect, log: &dyn LogSink) -> Result<NwkParse, Report> {
  let parse = tree_read_file(path, dialect)?;
  warn_duplicate_names(
    log,
    RunWarningKind::DuplicateNodeNames,
    &format!(
      "The tree '{}' gives the same name to more than one node:",
      path.display()
    ),
    "TreeTime keeps the name of the first such node and writes the others as <name>_2, <name>_3 and so on. Nodes with the same input name still receive the same data from the other inputs.",
    &parse.duplicate_names,
  );
  Ok(parse)
}

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct TreeArgs {
  /// Tree in Newick format
  ///
  /// A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
  #[cfg_attr(
    feature = "clap",
    clap(long, short = 't', value_hint = ValueHint::FilePath, help_heading = "Input data")
  )]
  #[schemars(extend("x-path" = "input"))]
  pub tree: Option<PathBuf>,
}

/// Segment trees and their alignments, for genomes whose segments have different trees.
#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct SegmentArgs {
  /// Trees of the genome segments, one Newick or Nexus file per segment, merged into one graph by node name
  ///
  /// Every node of every tree needs a name, unique within its tree; a name that occurs in several trees is one node
  /// of the graph. A node that a segment only passes through stays in that segment's tree as a node with one child.
  /// An edge shared by several trees gets the mean of their branch lengths. Use instead of `--tree`, together with
  /// `--segment-alignments`. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
  #[cfg_attr(
    feature = "clap",
    clap(
      long,
      num_args = 1..,
      value_delimiter = ',',
      value_hint = ValueHint::FilePath,
      hide_short_help = true,
      help_heading = "Input data"
    )
  )]
  #[schemars(extend("x-path" = "input"))]
  pub segment_trees: Vec<PathBuf>,

  /// Aligned FASTA input of the genome segments, one file per `--segment-trees` file, in the same order
  ///
  /// Use instead of `--alignment`. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
  #[cfg_attr(
    feature = "clap",
    clap(
      long,
      num_args = 1..,
      value_delimiter = ',',
      value_hint = ValueHint::FilePath,
      hide_short_help = true,
      help_heading = "Input data"
    )
  )]
  #[schemars(extend("x-path" = "input"))]
  pub segment_alignments: Vec<PathBuf>,
}

/// Newick dialect of the tree inputs, shared by every command that reads a tree.
#[derive(Debug, Clone, Copy, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct TreeDialectArgs {
  /// Newick dialect of the tree inputs (the tree and the reference topology)
  ///
  /// A dialect is a tree structure and an annotation convention, written as `<structure>,<annotations>`. The
  /// structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as `x#H1`) or rich (eNewick
  /// with support and probability fields). The annotations are plain (comments are text), beast (`[&key=value]`), nhx
  /// (`[&&NHX:key=value]`) or mrbayes (`[&B name value]`). For example, `enewick,beast` reads networks with BEAST
  /// annotations and `classic,nhx` reads NHX trees. TreeTime reads trees only, so a network input is an error.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, value_enum, default_value_t = TreeDialectArg::default(), help_heading = "Input data")
  )]
  #[deser(as = DisplayFromStr)]
  pub newick_dialect: TreeDialectArg,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Display, FromStr)]
pub struct TreeDialectArg(NewickDialect);

impl TreeDialectArg {
  pub fn dialect(self) -> NewickDialect {
    self.0
  }
}

impl Default for TreeDialectArg {
  fn default() -> Self {
    Self(TREE_DIALECT_DEFAULT)
  }
}

impl JsonSchema for TreeDialectArg {
  fn inline_schema() -> bool {
    true
  }

  fn schema_name() -> Cow<'static, str> {
    Cow::Borrowed("TreeDialectArg")
  }

  fn json_schema(_generator: &mut SchemaGenerator) -> Schema {
    json_schema!({
      "type": "string",
      "enum": TREE_DIALECT_NAMES.as_slice(),
    })
  }
}

#[cfg(feature = "clap")]
impl clap::ValueEnum for TreeDialectArg {
  fn value_variants<'a>() -> &'a [Self] {
    TREE_DIALECTS.as_slice()
  }

  fn to_possible_value(&self) -> Option<PossibleValue> {
    let idx = TREE_DIALECTS.iter().position(|dialect| dialect == self)?;
    TREE_DIALECT_NAMES
      .get(idx)
      .map(|name| PossibleValue::new(name.as_str()))
  }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TreeInput {
  Tree(PathBuf),
  Segments(Vec<PathBuf>),
}

impl TreeInput {
  pub fn tree(&self) -> Option<&Path> {
    match self {
      Self::Tree(tree) => Some(tree),
      Self::Segments(_) => None,
    }
  }

  pub fn is_segments(&self) -> bool {
    matches!(self, Self::Segments(_))
  }
}

pub(crate) fn resolve_tree_input(
  tree: Option<PathBuf>,
  alignment: Vec<PathBuf>,
  segments: SegmentArgs,
  missing_tree: impl FnOnce() -> Report,
) -> Result<(TreeInput, Vec<PathBuf>), Report> {
  let SegmentArgs {
    segment_trees,
    segment_alignments,
  } = segments;
  if segment_trees.is_empty() {
    if !segment_alignments.is_empty() {
      return Err(argument_error(
        &["segment_alignments", "segment_trees"],
        "--segment-alignments requires --segment-trees, one tree per alignment",
      ));
    }
    let tree = tree.ok_or_else(missing_tree)?;
    return Ok((TreeInput::Tree(tree), alignment));
  }
  if tree.is_some() {
    return Err(argument_error(
      &["tree", "segment_trees"],
      "--segment-trees replaces --tree: give one of them",
    ));
  }
  if !alignment.is_empty() {
    return Err(argument_error(
      &["alignment", "segment_alignments"],
      "--segment-trees takes its alignments from --segment-alignments, which replaces --alignment: give one of them",
    ));
  }
  if segment_alignments.len() != segment_trees.len() {
    return Err(argument_error(
      &["segment_trees", "segment_alignments"],
      format!(
        "--segment-trees and --segment-alignments pair one tree with one alignment, but they list {} trees and {} alignments",
        segment_trees.len(),
        segment_alignments.len()
      ),
    ));
  }
  Ok((TreeInput::Segments(segment_trees), segment_alignments))
}
