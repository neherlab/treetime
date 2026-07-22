#[cfg(feature = "clap")]
use clap::ValueHint;
use eyre::{Report, WrapErr};
use maplit::btreeset;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::fmt::Debug;
use std::path::{Path, PathBuf};
use strum::IntoEnumIterator;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_graph::topology_order::{TopologyOrderPreset, TopologyOrderSpec, TopologyOrderTargetAggregate};
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::NwkStyle;
use treetime_io::nwk::{EdgeFromNwk, NodeFromNwk, nwk_read_file};
use treetime_utils::{make_error, make_report};

/// Internal resolution and lookup key for every selectable output.
///
/// This is the lingua franca of the output system: per-command CLI selection enums convert into
/// it (`From<XxxOutputSelection>`), `resolve` maps it to concrete paths, and command code looks up
/// produced files by it. NWK annotation style is orthogonal and lives on `--output-nwk-style`, so
/// the tree variants here are style-agnostic (`Nwk`, `Nexus`), not per-style.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd, Serialize, Deserialize, strum_macros::EnumIter)]
pub enum OutputSelection {
  All,

  // Tree formats
  Nwk,
  Nexus,
  Auspice,
  Phyloxml,
  PhyloxmlJson,
  MatPb,
  MatJson,
  GraphJson,
  Dot,

  // Non-tree outputs
  AugurNodeData,
  Gtr,
  ClockModel,
  ConfidenceTsv,
  ConfidenceCsv,
  ReconstructedNucFasta,
  ReconstructedAaFasta,
  TraitsCsv,
  ClockCsv,
  Tracelog,
}

impl OutputSelection {
  pub fn is_tree(self) -> bool {
    matches!(
      self,
      Self::Nwk
        | Self::Nexus
        | Self::Auspice
        | Self::Phyloxml
        | Self::PhyloxmlJson
        | Self::MatPb
        | Self::MatJson
        | Self::GraphJson
        | Self::Dot
    )
  }

  /// Tree formats whose serialization is parameterized by NWK annotation style.
  pub fn is_styled_tree(self) -> bool {
    matches!(self, Self::Nwk | Self::Nexus)
  }

  pub fn is_meta(self) -> bool {
    matches!(self, Self::All)
  }

  pub fn extension(self) -> &'static str {
    match self {
      Self::All => "",
      Self::Nwk => ".nwk",
      Self::Nexus => ".nexus",
      Self::Auspice => ".auspice.json",
      Self::Phyloxml => ".phylo.xml",
      Self::PhyloxmlJson => ".phylo.json",
      Self::MatPb => ".mat.pb",
      Self::MatJson => ".mat.json",
      Self::GraphJson => ".graph.json",
      Self::Dot => ".dot",
      Self::AugurNodeData => ".augur-node-data.json",
      Self::Gtr => ".gtr.json",
      Self::ClockModel => ".clock-model.json",
      Self::ConfidenceTsv => ".confidence.tsv",
      Self::ConfidenceCsv => ".confidence.csv",
      Self::ReconstructedNucFasta => ".reconstructed-nuc.fasta",
      Self::ReconstructedAaFasta => ".reconstructed-aa.{cds}.fasta",
      Self::TraitsCsv => ".traits.csv",
      Self::ClockCsv => ".clock.csv",
      Self::Tracelog => ".tracelog.csv",
    }
  }

  /// Dispatch tag for the non-styled tree formats. Styled formats (`Nwk`, `Nexus`) carry a style
  /// and are converted via `styled_tree_write_kind`, so they return `None` here.
  pub fn to_tree_write_kind(self) -> Option<TreeWriteKind> {
    match self {
      Self::Auspice => Some(TreeWriteKind::Auspice),
      Self::Phyloxml => Some(TreeWriteKind::Phyloxml),
      Self::PhyloxmlJson => Some(TreeWriteKind::PhyloxmlJson),
      Self::MatPb => Some(TreeWriteKind::MatPb),
      Self::MatJson => Some(TreeWriteKind::MatJson),
      Self::GraphJson => Some(TreeWriteKind::GraphJson),
      Self::Dot => Some(TreeWriteKind::Dot),
      _ => None,
    }
  }

  pub fn flag_name(self) -> &'static str {
    match self {
      Self::All => "--output-selection=all",
      Self::Nwk => "--output-tree-nwk",
      Self::Nexus => "--output-tree-nexus",
      Self::Auspice => "--output-tree-auspice",
      Self::Phyloxml => "--output-tree-phyloxml",
      Self::PhyloxmlJson => "--output-tree-phyloxml-json",
      Self::MatPb => "--output-tree-mat-pb",
      Self::MatJson => "--output-tree-mat-json",
      Self::GraphJson => "--output-tree-graph-json",
      Self::Dot => "--output-tree-dot",
      Self::AugurNodeData => "--output-augur-node-data",
      Self::Gtr => "--output-gtr",
      Self::ClockModel => "--output-clock-model",
      Self::ConfidenceTsv => "--output-confidence-tsv",
      Self::ConfidenceCsv => "--output-confidence-csv",
      Self::ReconstructedNucFasta => "--output-reconstructed-nuc-fasta",
      Self::ReconstructedAaFasta => "--output-reconstructed-aa-fasta",
      Self::TraitsCsv => "--output-traits-csv",
      Self::ClockCsv => "--output-clock-csv",
      Self::Tracelog => "--output-tracelog",
    }
  }

  #[allow(clippy::option_option)]
  fn tree_field(self, args: &OutputCoreArgs) -> Option<Option<&Path>> {
    match self {
      Self::Nwk => Some(args.output_tree_nwk.as_deref()),
      Self::Nexus => Some(args.output_tree_nexus.as_deref()),
      Self::Auspice => Some(args.output_tree_auspice.as_deref()),
      Self::Phyloxml => Some(args.output_tree_phyloxml.as_deref()),
      Self::PhyloxmlJson => Some(args.output_tree_phyloxml_json.as_deref()),
      Self::MatPb => Some(args.output_tree_mat_pb.as_deref()),
      Self::MatJson => Some(args.output_tree_mat_json.as_deref()),
      Self::GraphJson => Some(args.output_tree_graph_json.as_deref()),
      Self::Dot => Some(args.output_tree_dot.as_deref()),
      _ => None,
    }
  }
}

impl std::fmt::Display for OutputSelection {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", self.flag_name())
  }
}

/// CLI-facing NWK/Nexus annotation style for `--output-nwk-style`.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
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

/// Secondary filename extension that distinguishes style-specific tree files when more than one
/// style is requested. Plain keeps the base name (no secondary extension).
fn nwk_style_secondary_ext(style: NwkStyle) -> &'static str {
  match style {
    NwkStyle::Plain => "",
    NwkStyle::Beast => ".annotated",
    NwkStyle::Nhx => ".nhx",
  }
}

/// Generates a per-command CLI selection enum plus its conversion to the internal `OutputSelection`.
///
/// Every command exposes the full tree-format surface (`Nwk`..`Dot`) plus `All`; the per-command
/// extras are the non-tree outputs that command supports. The CLI enum is the validation layer:
/// clap rejects values outside the command's variant set.
macro_rules! per_command_output_selection {
  ($name:ident { $($extra:ident),* $(,)? }) => {
    #[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
    #[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
    pub enum $name {
      All,
      Nwk,
      Nexus,
      Auspice,
      Phyloxml,
      PhyloxmlJson,
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
          $name::Phyloxml => Self::Phyloxml,
          $name::PhyloxmlJson => Self::PhyloxmlJson,
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
  ClockModel,
  ConfidenceTsv,
  Tracelog,
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

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum CommandKind {
  Ancestral,
  Timetree,
  Optimize,
  Mugration,
  Clock,
  Prune,
}

impl CommandKind {
  /// Full set selectable on this command.
  pub fn all_selectable(self) -> BTreeSet<OutputSelection> {
    &Self::available_tree_outputs() | &self.non_tree_outputs()
  }

  /// Outputs produced by `--output-all` without an explicit `--output-selection`.
  #[allow(clippy::enum_glob_use)]
  pub fn default_outputs(self) -> BTreeSet<OutputSelection> {
    use OutputSelection::*;
    let tree_defaults = if self == Self::Timetree {
      btreeset![Nwk, Nexus, Auspice]
    } else {
      btreeset![Nwk, Nexus]
    };
    let mut non_tree = self.non_tree_outputs();
    for non_default in [ReconstructedAaFasta, ConfidenceTsv, ConfidenceCsv, Tracelog] {
      non_tree.remove(&non_default);
    }
    &tree_defaults | &non_tree
  }

  #[allow(clippy::enum_glob_use)]
  fn available_tree_outputs() -> BTreeSet<OutputSelection> {
    use OutputSelection::*;
    btreeset![
      Nwk,
      Nexus,
      Auspice,
      Phyloxml,
      PhyloxmlJson,
      MatPb,
      MatJson,
      GraphJson,
      Dot
    ]
  }

  #[allow(clippy::enum_glob_use)]
  fn non_tree_outputs(self) -> BTreeSet<OutputSelection> {
    use OutputSelection::*;
    match self {
      Self::Ancestral => btreeset![AugurNodeData, Gtr, ReconstructedNucFasta, ReconstructedAaFasta],
      Self::Timetree => btreeset![AugurNodeData, Gtr, ClockModel, ConfidenceTsv, Tracelog],
      Self::Optimize => btreeset![AugurNodeData, Gtr],
      Self::Mugration => btreeset![AugurNodeData, Gtr, ConfidenceCsv, TraitsCsv],
      Self::Clock => btreeset![ClockModel, ClockCsv],
      Self::Prune => btreeset![Gtr],
    }
  }

  pub fn stem(self) -> &'static str {
    match self {
      Self::Ancestral => "ancestral",
      Self::Timetree => "timetree",
      Self::Optimize => "optimize",
      Self::Mugration => "mugration",
      Self::Clock => "clock",
      Self::Prune => "prune",
    }
  }

  /// Default NWK annotation styles when `--output-nwk-style` is not given. Plain for all commands;
  /// this is the extension point for per-command style defaults.
  pub fn default_nwk_styles(self) -> Vec<NwkStyle> {
    vec![NwkStyle::Plain]
  }
}

/// Insert a secondary filename extension before the final extension of a path.
/// `my.nwk` + `.annotated` -> `my.annotated.nwk`. Empty secondary leaves the path unchanged.
fn insert_secondary_ext(path: &Path, secondary: &str) -> PathBuf {
  if secondary.is_empty() {
    return path.to_path_buf();
  }
  if let (Some(stem), Some(ext)) = (path.file_stem(), path.extension()) {
    let mut name = OsString::from(stem);
    name.push(secondary);
    name.push(".");
    name.push(ext);
    path.with_file_name(name)
  } else {
    let mut name = OsString::from(path.as_os_str());
    name.push(secondary);
    PathBuf::from(name)
  }
}

fn styled_tree_write_kind(variant: OutputSelection, style: NwkStyle) -> TreeWriteKind {
  match variant {
    OutputSelection::Nwk => TreeWriteKind::nwk(style),
    OutputSelection::Nexus => TreeWriteKind::nexus(style),
    _ => unreachable!("styled_tree_write_kind called on non-styled variant"),
  }
}

/// Expand a single per-file NWK/Nexus override path across the selected styles.
///
/// A single style uses the override path verbatim; multiple styles insert each style's secondary
/// extension so the files do not collide.
fn expand_override_styles(path: &Path, styles: &[NwkStyle]) -> Vec<(NwkStyle, PathBuf)> {
  let multi = styles.len() > 1;
  styles
    .iter()
    .map(|&style| {
      let p = if multi {
        insert_secondary_ext(path, nwk_style_secondary_ext(style))
      } else {
        path.to_path_buf()
      };
      (style, p)
    })
    .collect()
}

/// Expand an `--output-all` NWK/Nexus output across the selected styles, deriving `{stem}{sec}{ext}`.
fn expand_outputall_styles(
  dir: &Path,
  stem: &str,
  variant: OutputSelection,
  styles: &[NwkStyle],
) -> Vec<(NwkStyle, PathBuf)> {
  let multi = styles.len() > 1;
  let primary = variant.extension();
  styles
    .iter()
    .map(|&style| {
      let secondary = if multi { nwk_style_secondary_ext(style) } else { "" };
      (style, dir.join(format!("{stem}{secondary}{primary}")))
    })
    .collect()
}

/// Three-tier output selection shared by every tree-writing command.
///
/// Tier 1: `--output-all` bulk directory with default file names.
/// Tier 2: `--output-selection` (a per-command field) restricts which outputs tier 1 produces.
/// Tier 3: Per-file `--output-tree-*` flags override or supplement tiers 1-2.
///
/// NWK annotation style (`--output-nwk-style`) is orthogonal and expands every NWK/Nexus output
/// across the selected styles. Topology ordering is a separate concern (`TopologyOrderArgs`) that
/// each command flattens independently.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
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

  /// Path to output PhyloXML tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_tree_phyloxml: Option<PathBuf>,

  /// Path to output PhyloXML-JSON tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  pub output_tree_phyloxml_json: Option<PathBuf>,

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

pub struct ResolvedOutputs {
  pub tree_outputs: BTreeMap<TreeWriteKind, PathBuf>,
  pub non_tree_outputs: BTreeMap<OutputSelection, PathBuf>,
}

impl OutputCoreArgs {
  /// Resolve the three-tier output configuration into concrete file paths.
  ///
  /// `selection` is the command's `--output-selection` already converted to `OutputSelection`.
  /// `non_tree_fields` carries the command's per-file non-tree flag values keyed by selection.
  ///
  /// Per-file flags are honored unconditionally and take precedence over `--output-all`. Every tree
  /// format is available to every command; runtime data prerequisites for non-tree outputs (e.g. a
  /// fitted GTR model) are checked by each command at write time, not here.
  pub fn resolve(
    &self,
    command: CommandKind,
    selection: &[OutputSelection],
    non_tree_fields: &[(OutputSelection, Option<&Path>)],
  ) -> Result<ResolvedOutputs, Report> {
    let stem = command.stem();
    let styles = self.effective_nwk_styles(command);

    let mut tree_outputs: BTreeMap<TreeWriteKind, PathBuf> = BTreeMap::new();
    let mut non_tree_outputs: BTreeMap<OutputSelection, PathBuf> = BTreeMap::new();

    // Tier 3a: per-file tree overrides, honored regardless of --output-all.
    let mut overridden_tree: BTreeSet<OutputSelection> = BTreeSet::new();
    for variant in OutputSelection::iter() {
      if !variant.is_tree() {
        continue;
      }
      let Some(Some(path)) = variant.tree_field(self) else {
        continue;
      };
      overridden_tree.insert(variant);
      if variant.is_styled_tree() {
        for (style, p) in expand_override_styles(path, &styles) {
          tree_outputs.insert(styled_tree_write_kind(variant, style), p);
        }
      } else {
        let kind = variant.to_tree_write_kind().expect("non-styled tree variant");
        tree_outputs.insert(kind, path.to_path_buf());
      }
    }

    // Tier 3b: per-file non-tree overrides.
    for &(sel, ref_path) in non_tree_fields {
      if let Some(path) = ref_path {
        non_tree_outputs.insert(sel, path.to_path_buf());
      }
    }

    // Tier 1+2: --output-all fills defaults or the explicit selection.
    if let Some(dir) = &self.output_all {
      let effective: BTreeSet<OutputSelection> = if selection.is_empty() {
        command.default_outputs()
      } else if selection.contains(&OutputSelection::All) {
        command.all_selectable()
      } else {
        selection.iter().copied().collect()
      };

      for &variant in &effective {
        if variant.is_meta() {
          continue;
        }
        if variant.is_tree() {
          if overridden_tree.contains(&variant) {
            continue;
          }
          if variant.is_styled_tree() {
            for (style, p) in expand_outputall_styles(dir, stem, variant, &styles) {
              tree_outputs.entry(styled_tree_write_kind(variant, style)).or_insert(p);
            }
          } else {
            let kind = variant.to_tree_write_kind().expect("non-styled tree variant");
            tree_outputs
              .entry(kind)
              .or_insert_with(|| dir.join(format!("{stem}{}", variant.extension())));
          }
        } else {
          non_tree_outputs
            .entry(variant)
            .or_insert_with(|| dir.join(format!("{stem}{}", variant.extension())));
        }
      }
    } else if !selection.is_empty() {
      return make_error!("--output-selection requires --output-all");
    }

    if tree_outputs.is_empty() && non_tree_outputs.is_empty() {
      return make_error!(
        "No output flags provided. At least one is required: \
         --output-all or one of the --output-tree-* / --output-* flags"
      );
    }

    ensure_unique_output_paths(&tree_outputs, &non_tree_outputs)?;

    Ok(ResolvedOutputs {
      tree_outputs,
      non_tree_outputs,
    })
  }

  /// Selected styles, de-duplicated and order-preserving, falling back to the command default.
  fn effective_nwk_styles(&self, command: CommandKind) -> Vec<NwkStyle> {
    if self.output_nwk_style.is_empty() {
      return command.default_nwk_styles();
    }
    let mut seen: BTreeSet<NwkStyle> = BTreeSet::new();
    let mut styles = Vec::new();
    for &style in &self.output_nwk_style {
      let style: NwkStyle = style.into();
      if seen.insert(style) {
        styles.push(style);
      }
    }
    styles
  }
}

fn ensure_unique_output_paths(
  tree_outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  non_tree_outputs: &BTreeMap<OutputSelection, PathBuf>,
) -> Result<(), Report> {
  let mut destinations: BTreeMap<&Path, String> = BTreeMap::new();
  for (kind, path) in tree_outputs {
    if let Some(previous) = destinations.insert(path, format!("{kind:?}")) {
      return make_error!(
        "Output destination '{}' is selected more than once ({previous} and {kind:?})",
        path.display()
      );
    }
  }
  for (selection, path) in non_tree_outputs {
    if let Some(previous) = destinations.insert(path, selection.flag_name().to_owned()) {
      return make_error!(
        "Output destination '{}' is selected more than once ({previous} and {})",
        path.display(),
        selection.flag_name()
      );
    }
  }
  Ok(())
}

#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct TopologyOrderArgs {
  /// Order tree topology before writing output files.
  #[cfg_attr(feature = "clap", clap(long, value_enum, help_heading = "Tree ordering"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order_target_source"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order_target_file"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order_target_aggregate"))]
  pub ladderize: Option<LadderizeArg>,

  /// Canonical topology ordering preset.
  #[cfg_attr(feature = "clap", clap(long, value_enum, help_heading = "Tree ordering"))]
  pub topology_order: Option<TopologyOrderArg>,

  /// Source for target-order topology sorting.
  #[cfg_attr(feature = "clap", clap(long, value_enum, help_heading = "Tree ordering"))]
  pub topology_order_target_source: Option<TopologyOrderTargetSourceArg>,

  /// File used by list or reference-topology target-order sources.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Tree ordering"))]
  pub topology_order_target_file: Option<PathBuf>,

  /// Aggregate used to map a subtree to a target-order position.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = TopologyOrderTargetAggregateArg::default(), help_heading = "Tree ordering"))]
  #[default(TopologyOrderTargetAggregateArg::Mean)]
  pub topology_order_target_aggregate: TopologyOrderTargetAggregateArg,
}

impl TopologyOrderArgs {
  pub fn resolve_topology_order<N, E, D>(
    &self,
    graph: &Graph<N, E, D>,
    input_order: Option<Vec<String>>,
  ) -> Result<TopologyOrderSpec, Report>
  where
    N: GraphNode + Named,
    E: GraphEdge,
    D: Sync + Send,
  {
    self.validate()?;

    let preset = match (self.ladderize, self.topology_order) {
      (Some(LadderizeArg::None), None) => TopologyOrderPreset::Keep,
      (Some(LadderizeArg::Ascending), None) => TopologyOrderPreset::DescendantCount,
      (Some(LadderizeArg::Descending), None) => TopologyOrderPreset::DescendantCountReverse,
      (None, Some(topology_order)) => topology_order.into(),
      (None, None) => TopologyOrderPreset::DescendantCount,
      (Some(_), Some(_)) => {
        return make_error!("--ladderize cannot be combined with --topology-order");
      },
    };

    let target_order = if preset.is_target_order() {
      self.target_order(graph, input_order)?
    } else {
      vec![]
    };

    Ok(TopologyOrderSpec {
      preset,
      target_order,
      target_aggregate: self.topology_order_target_aggregate.into(),
    })
  }

  fn validate(&self) -> Result<(), Report> {
    if self.ladderize.is_some()
      && (self.topology_order.is_some()
        || self.topology_order_target_source.is_some()
        || self.topology_order_target_file.is_some()
        || self.topology_order_target_aggregate != TopologyOrderTargetAggregateArg::Mean)
    {
      return make_error!("--ladderize cannot be combined with --topology-order* options");
    }

    let target_fields_present = self.topology_order_target_source.is_some()
      || self.topology_order_target_file.is_some()
      || self.topology_order_target_aggregate != TopologyOrderTargetAggregateArg::Mean;
    let target_mode = self
      .topology_order
      .is_some_and(|order| order.into_preset().is_target_order());

    if target_fields_present && !target_mode {
      return make_error!(
        "--topology-order-target-* options require --topology-order=target-order or target-order-reverse"
      );
    }

    if target_mode
      && matches!(
        self.topology_order_target_source,
        Some(TopologyOrderTargetSourceArg::ReferenceTopology | TopologyOrderTargetSourceArg::List)
      )
      && self.topology_order_target_file.is_none()
    {
      return make_error!("--topology-order-target-file is required for reference-topology and list target sources");
    }

    Ok(())
  }

  fn target_order<N, E, D>(
    &self,
    graph: &Graph<N, E, D>,
    input_order: Option<Vec<String>>,
  ) -> Result<Vec<String>, Report>
  where
    N: GraphNode + Named,
    E: GraphEdge,
    D: Sync + Send,
  {
    match self
      .topology_order_target_source
      .unwrap_or(TopologyOrderTargetSourceArg::Input)
    {
      TopologyOrderTargetSourceArg::Input => input_order.map_or_else(|| leaf_order(graph), Ok),
      TopologyOrderTargetSourceArg::ReferenceTopology => {
        let path = self
          .topology_order_target_file
          .as_ref()
          .ok_or_else(|| make_report!("--topology-order-target-file is required for reference-topology"))?;
        let graph =
          nwk_read_file::<OrderNode, OrderEdge, ()>(path).wrap_err("When reading target reference topology")?;
        leaf_order(&graph)
      },
      TopologyOrderTargetSourceArg::List => {
        let path = self
          .topology_order_target_file
          .as_ref()
          .ok_or_else(|| make_report!("--topology-order-target-file is required for list"))?;
        let contents = std::fs::read_to_string(path)
          .wrap_err_with(|| format!("When reading topology order target list '{}'", path.display()))?;
        Ok(
          contents
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .map(str::to_owned)
            .collect(),
        )
      },
    }
  }
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
pub enum DivergenceUnits {
  #[default]
  MutationsPerSite,
  Mutations,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
pub enum LadderizeArg {
  None,
  Ascending,
  Descending,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
pub enum TopologyOrderArg {
  Keep,
  #[cfg_attr(feature = "clap", value(alias = "ladderize"))]
  DescendantCount,
  #[cfg_attr(feature = "clap", value(alias = "ladderize-reverse"))]
  DescendantCountReverse,
  Height,
  HeightReverse,
  Divergence,
  DivergenceReverse,
  Label,
  #[cfg_attr(feature = "clap", value(alias = "alphabetical-reverse"))]
  LabelReverse,
  TargetOrder,
  TargetOrderReverse,
}

impl TopologyOrderArg {
  fn into_preset(self) -> TopologyOrderPreset {
    self.into()
  }
}

impl From<TopologyOrderArg> for TopologyOrderPreset {
  fn from(value: TopologyOrderArg) -> Self {
    match value {
      TopologyOrderArg::Keep => Self::Keep,
      TopologyOrderArg::DescendantCount => Self::DescendantCount,
      TopologyOrderArg::DescendantCountReverse => Self::DescendantCountReverse,
      TopologyOrderArg::Height => Self::Height,
      TopologyOrderArg::HeightReverse => Self::HeightReverse,
      TopologyOrderArg::Divergence => Self::Divergence,
      TopologyOrderArg::DivergenceReverse => Self::DivergenceReverse,
      TopologyOrderArg::Label => Self::Label,
      TopologyOrderArg::LabelReverse => Self::LabelReverse,
      TopologyOrderArg::TargetOrder => Self::TargetOrder,
      TopologyOrderArg::TargetOrderReverse => Self::TargetOrderReverse,
    }
  }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
pub enum TopologyOrderTargetSourceArg {
  Input,
  ReferenceTopology,
  List,
}

impl std::fmt::Display for TopologyOrderTargetSourceArg {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    match self {
      Self::Input => write!(f, "input"),
      Self::ReferenceTopology => write!(f, "reference-topology"),
      Self::List => write!(f, "list"),
    }
  }
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
pub enum TopologyOrderTargetAggregateArg {
  #[default]
  Mean,
  Median,
}

impl From<TopologyOrderTargetAggregateArg> for TopologyOrderTargetAggregate {
  fn from(value: TopologyOrderTargetAggregateArg) -> Self {
    match value {
      TopologyOrderTargetAggregateArg::Mean => Self::Mean,
      TopologyOrderTargetAggregateArg::Median => Self::Median,
    }
  }
}

fn leaf_order<N, E, D>(graph: &Graph<N, E, D>) -> Result<Vec<String>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  graph
    .get_leaves()
    .into_iter()
    .map(|leaf| {
      let leaf = leaf.read_arc();
      leaf
        .payload()
        .read_arc()
        .name()
        .map(|name| name.as_ref().to_owned())
        .ok_or_else(|| make_report!("Leaf node {} has no name", leaf.key()))
    })
    .collect()
}

#[derive(Clone, Debug, Default)]
struct OrderNode {
  name: Option<String>,
}

impl GraphNode for OrderNode {}

impl Named for OrderNode {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|name| name.as_ref().to_owned());
  }
}

impl NodeFromNwk for OrderNode {
  fn from_nwk(
    name: Option<impl AsRef<str>>,
    _confidence: Option<f64>,
    _: &BTreeMap<String, String>,
  ) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|name| name.as_ref().to_owned()),
    })
  }
}

#[derive(Clone, Debug, Default)]
struct OrderEdge {
  branch_length: Option<f64>,
}

impl GraphEdge for OrderEdge {}

impl HasBranchLength for OrderEdge {
  fn branch_length(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_branch_length(&mut self, branch_length: Option<f64>) {
    self.branch_length = branch_length;
  }
}

impl EdgeFromNwk for OrderEdge {
  fn from_nwk(weight: Option<f64>) -> Result<Self, Report> {
    Ok(Self { branch_length: weight })
  }
}
