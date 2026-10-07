use crate::commands::shared::input_warnings::warn_duplicate_names;
#[cfg(feature = "clap")]
use clap::builder::PossibleValue;
use derive_more::{Display, FromStr};
use deser::adapters::DisplayFromStr;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::{JsonSchema, Schema, SchemaGenerator, json_schema};
use smart_default::SmartDefault;
use std::borrow::Cow;
use std::path::Path;
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
    "Nodes with the same name receive the same data from the other inputs, and augur node data keeps one entry per name.",
    &parse.duplicate_names,
  );
  Ok(parse)
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
  /// A dialect is a tree structure and an annotation convention, written as <structure>,<annotations>. The structure
  /// is classic (plain labels and branch lengths), enewick (hybrid nodes such as x#H1) or rich (eNewick with support
  /// and probability fields). The annotations are plain (comments are text), beast ([&key=value]), nhx
  /// ([&&NHX:key=value]) or mrbayes ([&B name value]). For example, enewick,beast reads networks with BEAST annotations
  /// and classic,nhx reads NHX trees. TreeTime reads trees only, so a network input is an error.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_enum, default_value_t = TreeDialectArg::default(), help_heading = "Input data")
  )]
  #[deser(as = DisplayFromStr)]
  pub tree_dialect: TreeDialectArg,
}

impl TreeDialectArgs {
  pub fn dialect(&self) -> NewickDialect {
    self.tree_dialect.dialect()
  }
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
