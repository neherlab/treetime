//! Output selection and path planning, shared by every application adapter.
//!
//! Turns a requested set of outputs plus a base destination into the concrete file paths and format
//! dispatch tags the writers consume. This tier is pure computation: it takes already-parsed values
//! (an [`OutputPlanRequest`]), opens no file, reads no input, and calls no format parser. An adapter
//! parses its own flags into the request, chooses the base destination, and opens the planned paths.

use eyre::Report;
use maplit::btreeset;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::path::{Path, PathBuf};
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::NwkStyle;
use treetime_utils::make_error;

/// Internal resolution and lookup key for every selectable output.
///
/// This is the lingua franca of the output system: per-command adapter selection enums convert into
/// it (`From<XxxOutputSelection>`), [`plan`] maps it to concrete paths, and command code looks up
/// produced files by it. NWK annotation style is orthogonal and lives on the adapter's
/// `--output-nwk-style` flag, so the tree variants here are style-agnostic (`Nwk`, `Nexus`), not
/// per-style.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum OutputSelection {
  All,

  // Tree formats
  Nwk,
  Nexus,
  Auspice,
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
  CoalescentTsv,
  CoalescentCsv,
  CoalescentJson,
}

impl OutputSelection {
  pub fn is_tree(self) -> bool {
    matches!(
      self,
      Self::Nwk | Self::Nexus | Self::Auspice | Self::MatPb | Self::MatJson | Self::GraphJson | Self::Dot
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
      Self::CoalescentTsv => ".coalescent.tsv",
      Self::CoalescentCsv => ".coalescent.csv",
      Self::CoalescentJson => ".coalescent.json",
    }
  }

  /// Dispatch tag for the non-styled tree formats. Styled formats (`Nwk`, `Nexus`) carry a style
  /// and are converted via `styled_tree_write_kind`, so they return `None` here.
  pub fn to_tree_write_kind(self) -> Option<TreeWriteKind> {
    match self {
      Self::Auspice => Some(TreeWriteKind::Auspice),
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
      Self::CoalescentTsv => "--output-coalescent-tsv",
      Self::CoalescentCsv => "--output-coalescent-csv",
      Self::CoalescentJson => "--output-coalescent-json",
    }
  }
}

impl std::fmt::Display for OutputSelection {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", self.flag_name())
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

/// The set of operations, each with its own selectable outputs and default file names.
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
    for non_default in [
      ReconstructedAaFasta,
      ConfidenceTsv,
      ConfidenceCsv,
      Tracelog,
      CoalescentCsv,
      CoalescentJson,
    ] {
      non_tree.remove(&non_default);
    }
    &tree_defaults | &non_tree
  }

  #[allow(clippy::enum_glob_use)]
  fn available_tree_outputs() -> BTreeSet<OutputSelection> {
    use OutputSelection::*;
    btreeset![Nwk, Nexus, Auspice, MatPb, MatJson, GraphJson, Dot]
  }

  #[allow(clippy::enum_glob_use)]
  fn non_tree_outputs(self) -> BTreeSet<OutputSelection> {
    use OutputSelection::*;
    match self {
      Self::Ancestral => btreeset![AugurNodeData, Gtr, ReconstructedNucFasta, ReconstructedAaFasta],
      Self::Timetree => btreeset![
        AugurNodeData,
        Gtr,
        ReconstructedNucFasta,
        ClockModel,
        ConfidenceTsv,
        Tracelog,
        CoalescentTsv,
        CoalescentCsv,
        CoalescentJson
      ],
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

  /// Default NWK annotation styles when the adapter passes no style. Plain for all commands; this is
  /// the extension point for per-command style defaults.
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

/// A requested set of outputs plus the base destination, ready for [`plan`] to turn into paths.
///
/// The adapter parses its flags into this value: it picks the base directory (`output_all`), the
/// requested annotation styles (empty means "use the command default"), the restricting selection,
/// and the per-file tree and non-tree destination overrides. Planning reads only these values and
/// performs no I/O.
pub struct OutputPlanRequest {
  /// Which command's output taxonomy and default file names apply.
  pub command: CommandKind,

  /// Base directory for `--output-all` bulk output, if the adapter requested one.
  pub output_all: Option<PathBuf>,

  /// Requested NWK/Nexus annotation styles. Empty means the command default.
  pub nwk_styles: Vec<NwkStyle>,

  /// Restricting selection for `--output-all`. Empty means the command's default set; `All` means
  /// the complete selectable set.
  pub selection: Vec<OutputSelection>,

  /// Per-file tree destination overrides, honored unconditionally and above `--output-all`.
  pub tree_overrides: BTreeMap<OutputSelection, PathBuf>,

  /// Per-file non-tree destination overrides, honored unconditionally and above `--output-all`.
  pub non_tree_overrides: BTreeMap<OutputSelection, PathBuf>,
}

/// Concrete output destinations, grouped by writer dispatch (tree formats) and selection (the rest).
pub struct ResolvedOutputs {
  pub tree_outputs: BTreeMap<TreeWriteKind, PathBuf>,
  pub non_tree_outputs: BTreeMap<OutputSelection, PathBuf>,
}

impl ResolvedOutputs {
  /// Group every produced file under the `OutputSelection` it satisfies.
  ///
  /// A styled tree format (`nwk`, `nexus`) expands to one path per requested annotation style, so a
  /// selection can map to several files. The pipeline uses this to resolve `{{ steps.x.outputs.<sel>
  /// }}` chaining and to reject a reference whose selection is ambiguous (more than one file). Paths
  /// within a selection are sorted for deterministic diagnostics.
  pub fn paths_by_selection(&self) -> BTreeMap<OutputSelection, Vec<PathBuf>> {
    let mut by_selection: BTreeMap<OutputSelection, Vec<PathBuf>> = BTreeMap::new();
    for (kind, path) in &self.tree_outputs {
      by_selection
        .entry(tree_write_kind_selection(kind))
        .or_default()
        .push(path.clone());
    }
    for (selection, path) in &self.non_tree_outputs {
      by_selection.entry(*selection).or_default().push(path.clone());
    }
    for paths in by_selection.values_mut() {
      paths.sort();
    }
    by_selection
  }
}

/// Invert the tree write dispatch tag back to the style-agnostic selection it was produced for.
fn tree_write_kind_selection(kind: &TreeWriteKind) -> OutputSelection {
  match kind {
    TreeWriteKind::Nwk(_) => OutputSelection::Nwk,
    TreeWriteKind::Nexus(_) => OutputSelection::Nexus,
    TreeWriteKind::Auspice => OutputSelection::Auspice,
    TreeWriteKind::MatPb => OutputSelection::MatPb,
    TreeWriteKind::MatJson => OutputSelection::MatJson,
    TreeWriteKind::GraphJson => OutputSelection::GraphJson,
    TreeWriteKind::Dot => OutputSelection::Dot,
  }
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
/// Resolve the three-tier output request into concrete file paths.
///
/// Tier 1: `output_all` bulk directory with default file names.
/// Tier 2: `selection` restricts which outputs tier 1 produces.
/// Tier 3: per-file `tree_overrides`/`non_tree_overrides` override or supplement tiers 1-2.
///
/// Per-file overrides are honored unconditionally and take precedence over `output_all`. Every tree
/// format is available to every command; runtime data prerequisites for non-tree outputs (e.g. a
/// fitted GTR model) are checked by each command at write time, not here.
///
/// NWK annotation style expands every NWK/Nexus output across the selected styles. Topology ordering
/// is a separate adapter concern applied to the graph before the writers run, not part of planning.
pub fn plan(request: &OutputPlanRequest) -> Result<ResolvedOutputs, Report> {
  let command = request.command;
  let stem = command.stem();
  let styles = effective_nwk_styles(command, &request.nwk_styles);

  let mut tree_outputs: BTreeMap<TreeWriteKind, PathBuf> = BTreeMap::new();
  let mut non_tree_outputs: BTreeMap<OutputSelection, PathBuf> = BTreeMap::new();

  // Tier 3a: per-file tree overrides, honored regardless of output_all.
  let mut overridden_tree: BTreeSet<OutputSelection> = BTreeSet::new();
  for (&variant, path) in &request.tree_overrides {
    overridden_tree.insert(variant);
    if variant.is_styled_tree() {
      for (style, p) in expand_override_styles(path, &styles) {
        tree_outputs.insert(styled_tree_write_kind(variant, style), p);
      }
    } else {
      let kind = variant.to_tree_write_kind().expect("non-styled tree variant");
      tree_outputs.insert(kind, path.clone());
    }
  }

  // Tier 3b: per-file non-tree overrides.
  for (&sel, path) in &request.non_tree_overrides {
    non_tree_outputs.insert(sel, path.clone());
  }

  // Tier 1+2: output_all fills defaults or the explicit selection.
  if let Some(dir) = &request.output_all {
    let effective: BTreeSet<OutputSelection> = if request.selection.is_empty() {
      command.default_outputs()
    } else if request.selection.contains(&OutputSelection::All) {
      command.all_selectable()
    } else {
      request.selection.iter().copied().collect()
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
  } else if !request.selection.is_empty() {
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
fn effective_nwk_styles(command: CommandKind, requested: &[NwkStyle]) -> Vec<NwkStyle> {
  if requested.is_empty() {
    return command.default_nwk_styles();
  }
  let mut seen: BTreeSet<NwkStyle> = BTreeSet::new();
  let mut styles = Vec::new();
  for &style in requested {
    if seen.insert(style) {
      styles.push(style);
    }
  }
  styles
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
