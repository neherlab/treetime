//! Command-line output arguments and their conversion into an output plan.
//!
//! Owns the CLI-facing output surface: the `--output-all`/`--output-selection`/`--output-tree-*`
//! flags, the per-command selection enums clap validates against, the `--output-nwk-style` flag, and
//! the `--divergence-units` flag. It parses these flags into an [`OutputPlanRequest`] and calls the
//! client-agnostic planner in `app_output::output_plan`, which computes the concrete paths. Opening
//! the planned files is each command's job.

use app_output::output_plan::{self, CommandKind, OutputPlanRequest, OutputSelection, ResolvedOutputs};
#[cfg(feature = "clap")]
use clap::ValueHint;
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime_io::nwk::NwkStyle;

/// CLI-facing NWK/Nexus annotation style for `--output-nwk-style`.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
pub enum NwkStyleArg {
  Plain,
  Beast,
  Nhx,
}

impl From<NwkStyleArg> for NwkStyle {
  fn from(value: NwkStyleArg) -> Self {
    match value {
      NwkStyleArg::Plain => Self::Plain,
      NwkStyleArg::Beast => Self::Beast,
      NwkStyleArg::Nhx => Self::Nhx,
    }
  }
}

/// Generates a per-command CLI selection enum plus its conversion to the internal `OutputSelection`.
///
/// Every command exposes the full tree-format surface (`Nwk`..`Dot`) plus `All`; the per-command
/// extras are the non-tree outputs that command supports. The CLI enum is the validation layer:
/// clap rejects values outside the command's variant set.
macro_rules! per_command_output_selection {
  ($name:ident { $($extra:ident),* $(,)? }) => {
    #[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize, schemars::JsonSchema)]
    #[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
    pub enum $name {
      All,
      Nwk,
      Nexus,
      Auspice,
      MatPb,
      MatJson,
      GraphJson,
      Dot,
      $($extra),*
    }

    impl From<$name> for OutputSelection {
      fn from(value: $name) -> Self {
        match value {
          $name::All => Self::All,
          $name::Nwk => Self::Nwk,
          $name::Nexus => Self::Nexus,
          $name::Auspice => Self::Auspice,
          $name::MatPb => Self::MatPb,
          $name::MatJson => Self::MatJson,
          $name::GraphJson => Self::GraphJson,
          $name::Dot => Self::Dot,
          $($name::$extra => Self::$extra),*
        }
      }
    }
  };
}

per_command_output_selection!(AncestralOutputSelection {
  AugurNodeData,
  Gtr,
  ReconstructedNucFasta,
  ReconstructedAaFasta,
});
per_command_output_selection!(TimetreeOutputSelection {
  AugurNodeData,
  Gtr,
  ReconstructedNucFasta,
  ClockModel,
  ConfidenceTsv,
  Tracelog,
  CoalescentTsv,
  CoalescentCsv,
  CoalescentJson,
});
per_command_output_selection!(ClockOutputSelection { ClockModel, ClockCsv });
per_command_output_selection!(MugrationOutputSelection {
  AugurNodeData,
  Gtr,
  ConfidenceCsv,
  TraitsCsv,
});
per_command_output_selection!(OptimizeOutputSelection { AugurNodeData, Gtr });
per_command_output_selection!(PruneOutputSelection { Gtr });

/// Three-tier output selection shared by every tree-writing command.
///
/// Tier 1: `--output-all` bulk directory with default file names.
/// Tier 2: `--output-selection` (a per-command field) restricts which outputs tier 1 produces.
/// Tier 3: Per-file `--output-tree-*` flags override or supplement tiers 1-2.
///
/// NWK annotation style (`--output-nwk-style`) is orthogonal and expands every NWK/Nexus output
/// across the selected styles. Topology ordering is a separate concern (`TopologyOrderArgs`) that
/// each command flattens independently.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct OutputCoreArgs {
  /// Write all default output files into this directory.
  ///
  /// Produces the default set of tree and non-tree outputs for the command, using
  /// `<dir>/<command>.<ext>` paths. Combine with `--output-selection` to restrict which
  /// outputs are written.
  ///
  /// Per-file flags (`--output-tree-nwk`, `--output-augur-node-data`, etc.) override or
  /// supplement the files produced by `--output-all`.
  #[cfg_attr(feature = "clap", clap(long, short = 'O', value_hint = ValueHint::DirPath, help_heading = "Output"))]
  pub output_all: Option<PathBuf>,

  /// NWK/Nexus annotation styles to write (comma-separated): `plain`, `beast`, `nhx`.
  ///
  /// Applies to every NWK and Nexus output. With more than one style, files are distinguished by a
  /// secondary extension (`.annotated` for beast, `.nhx` for nhx). Default: `plain`.
  #[cfg_attr(feature = "clap", clap(long, value_delimiter = ',', help_heading = "Output"))]
  pub output_nwk_style: Vec<NwkStyleArg>,

  /// Path to output Newick tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`. With
  /// multiple `--output-nwk-style` values, a secondary extension is inserted per style.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_tree_nwk: Option<PathBuf>,

  /// Path to output Nexus tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`. With
  /// multiple `--output-nwk-style` values, a secondary extension is inserted per style.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_tree_nexus: Option<PathBuf>,

  /// Path to output Auspice v2 JSON tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_tree_auspice: Option<PathBuf>,

  /// Path to output UShER MAT protobuf tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_tree_mat_pb: Option<PathBuf>,

  /// Path to output UShER MAT JSON tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_tree_mat_json: Option<PathBuf>,

  /// Path to output internal graph JSON tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_tree_graph_json: Option<PathBuf>,

  /// Path to output Graphviz DOT tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_tree_dot: Option<PathBuf>,
}

impl OutputCoreArgs {
  /// Resolve the three-tier output configuration into concrete file paths.
  ///
  /// `selection` is the command's `--output-selection` already converted to `OutputSelection`.
  /// `non_tree_fields` carries the command's per-file non-tree flag values keyed by selection.
  ///
  /// Parses the CLI flags into an [`OutputPlanRequest`] and delegates the path computation to the
  /// client-agnostic planner in `app_output::output_plan`.
  pub fn resolve(
    &self,
    command: CommandKind,
    selection: &[OutputSelection],
    non_tree_fields: &[(OutputSelection, Option<&Path>)],
  ) -> Result<ResolvedOutputs, Report> {
    let request = OutputPlanRequest {
      command,
      output_all: self.output_all.clone(),
      nwk_styles: self.output_nwk_style.iter().copied().map(NwkStyle::from).collect(),
      selection: selection.to_vec(),
      tree_overrides: self.tree_overrides(),
      non_tree_overrides: non_tree_fields
        .iter()
        .filter_map(|&(sel, path)| path.map(|path| (sel, path.to_path_buf())))
        .collect(),
    };
    output_plan::plan(&request)
  }

  /// Per-file tree destination overrides keyed by selection, one entry per `--output-tree-*` flag set.
  fn tree_overrides(&self) -> BTreeMap<OutputSelection, PathBuf> {
    [
      (OutputSelection::Nwk, self.output_tree_nwk.as_deref()),
      (OutputSelection::Nexus, self.output_tree_nexus.as_deref()),
      (OutputSelection::Auspice, self.output_tree_auspice.as_deref()),
      (OutputSelection::MatPb, self.output_tree_mat_pb.as_deref()),
      (OutputSelection::MatJson, self.output_tree_mat_json.as_deref()),
      (OutputSelection::GraphJson, self.output_tree_graph_json.as_deref()),
      (OutputSelection::Dot, self.output_tree_dot.as_deref()),
    ]
    .into_iter()
    .filter_map(|(variant, path)| path.map(|path| (variant, path.to_path_buf())))
    .collect()
  }
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
pub enum DivergenceUnits {
  #[default]
  MutationsPerSite,
  Mutations,
}
