use eyre::Report;
use maplit::btreeset;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::path::{Path, PathBuf};
use strum_macros::{AsRefStr, EnumIter, EnumString};
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::NwkStyle;
use treetime_utils::make_error;

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
pub fn plan(request: &OutputPlanRequest) -> Result<ResolvedOutputs, Report> {
  let command = request.command;
  let stem = command.stem();
  let styles = effective_nwk_styles(command, &request.nwk_styles);

  let mut tree_outputs: BTreeMap<TreeWriteKind, PathBuf> = BTreeMap::new();
  let mut non_tree_outputs: BTreeMap<OutputSelection, PathBuf> = BTreeMap::new();

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

  for (&sel, path) in &request.non_tree_overrides {
    non_tree_outputs.insert(sel, path.clone());
  }

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

fn styled_tree_write_kind(variant: OutputSelection, style: NwkStyle) -> TreeWriteKind {
  match variant {
    OutputSelection::Nwk => TreeWriteKind::nwk(style),
    OutputSelection::Nexus => TreeWriteKind::nexus(style),
    _ => unreachable!("styled_tree_write_kind called on non-styled variant"),
  }
}

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

fn nwk_style_secondary_ext(style: NwkStyle) -> &'static str {
  match style {
    NwkStyle::Plain => "",
    NwkStyle::Beast => ".annotated",
    NwkStyle::Nhx => ".nhx",
  }
}

pub struct OutputPlanRequest {
  pub command: CommandKind,

  pub output_all: Option<PathBuf>,

  pub nwk_styles: Vec<NwkStyle>,

  pub selection: Vec<OutputSelection>,

  pub tree_overrides: BTreeMap<OutputSelection, PathBuf>,

  pub non_tree_overrides: BTreeMap<OutputSelection, PathBuf>,
}

pub struct ResolvedOutputs {
  pub tree_outputs: BTreeMap<TreeWriteKind, PathBuf>,
  pub non_tree_outputs: BTreeMap<OutputSelection, PathBuf>,
}

impl ResolvedOutputs {
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
  pub fn all_selectable(self) -> BTreeSet<OutputSelection> {
    &Self::available_tree_outputs() | &self.non_tree_outputs()
  }

  fn default_outputs(self) -> BTreeSet<OutputSelection> {
    let tree_defaults = if self == Self::Timetree {
      btreeset![OutputSelection::Nwk, OutputSelection::Nexus, OutputSelection::Auspice]
    } else {
      btreeset![OutputSelection::Nwk, OutputSelection::Nexus]
    };
    let mut non_tree = self.non_tree_outputs();
    for non_default in [
      OutputSelection::ReconstructedAaFasta,
      OutputSelection::ConfidenceTsv,
      OutputSelection::ConfidenceCsv,
      OutputSelection::Tracelog,
      OutputSelection::CoalescentCsv,
      OutputSelection::CoalescentJson,
    ] {
      non_tree.remove(&non_default);
    }
    &tree_defaults | &non_tree
  }

  fn available_tree_outputs() -> BTreeSet<OutputSelection> {
    btreeset![
      OutputSelection::Nwk,
      OutputSelection::Nexus,
      OutputSelection::Auspice,
      OutputSelection::MatPb,
      OutputSelection::MatJson,
      OutputSelection::GraphJson,
      OutputSelection::Dot
    ]
  }

  fn non_tree_outputs(self) -> BTreeSet<OutputSelection> {
    match self {
      Self::Ancestral => btreeset![
        OutputSelection::AugurNodeData,
        OutputSelection::Gtr,
        OutputSelection::ReconstructedNucFasta,
        OutputSelection::ReconstructedAaFasta
      ],
      Self::Timetree => btreeset![
        OutputSelection::AugurNodeData,
        OutputSelection::Gtr,
        OutputSelection::ReconstructedNucFasta,
        OutputSelection::ClockModel,
        OutputSelection::ConfidenceTsv,
        OutputSelection::Tracelog,
        OutputSelection::CoalescentTsv,
        OutputSelection::CoalescentCsv,
        OutputSelection::CoalescentJson
      ],
      Self::Optimize => btreeset![OutputSelection::AugurNodeData, OutputSelection::Gtr],
      Self::Mugration => btreeset![
        OutputSelection::AugurNodeData,
        OutputSelection::Gtr,
        OutputSelection::ConfidenceCsv,
        OutputSelection::TraitsCsv
      ],
      Self::Clock => btreeset![OutputSelection::ClockModel, OutputSelection::ClockCsv],
      Self::Prune => btreeset![OutputSelection::Gtr],
    }
  }

  fn stem(self) -> &'static str {
    match self {
      Self::Ancestral => "ancestral",
      Self::Timetree => "timetree",
      Self::Optimize => "optimize",
      Self::Mugration => "mugration",
      Self::Clock => "clock",
      Self::Prune => "prune",
    }
  }

  fn default_nwk_styles(self) -> Vec<NwkStyle> {
    vec![NwkStyle::Plain]
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

/// Canonical lookup key for selectable outputs. Command adapters convert their
/// selection enums into this type, and [`plan`] resolves each key to a path.
/// Tree variants do not encode the separately selected Newick style.
#[derive(
  Copy,
  Clone,
  Debug,
  Eq,
  PartialEq,
  Hash,
  Ord,
  PartialOrd,
  Serialize,
  Deserialize,
  JsonSchema,
  AsRefStr,
  EnumString,
  EnumIter,
)]
#[serde(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum OutputSelection {
  All,

  Nwk,
  Nexus,
  Auspice,
  MatPb,
  MatJson,
  GraphJson,
  Dot,

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
  fn is_tree(self) -> bool {
    matches!(
      self,
      Self::Nwk | Self::Nexus | Self::Auspice | Self::MatPb | Self::MatJson | Self::GraphJson | Self::Dot
    )
  }

  fn is_styled_tree(self) -> bool {
    matches!(self, Self::Nwk | Self::Nexus)
  }

  fn is_meta(self) -> bool {
    matches!(self, Self::All)
  }

  fn extension(self) -> &'static str {
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

  fn to_tree_write_kind(self) -> Option<TreeWriteKind> {
    match self {
      Self::Auspice => Some(TreeWriteKind::Auspice),
      Self::MatPb => Some(TreeWriteKind::MatPb),
      Self::MatJson => Some(TreeWriteKind::MatJson),
      Self::GraphJson => Some(TreeWriteKind::GraphJson),
      Self::Dot => Some(TreeWriteKind::Dot),
      _ => None,
    }
  }

  fn flag_name(self) -> &'static str {
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
