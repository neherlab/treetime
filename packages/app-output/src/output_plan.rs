use derive_more::{Display, Error};
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use maplit::btreeset;
use schemars::JsonSchema;
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::num::NonZeroUsize;
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf, is_separator};
use strum_macros::{AsRefStr, EnumIter, EnumString};
use treetime::progress::LogSink;
use treetime::progress_debug;
use treetime_io::csv::TableFormat;
use treetime_io::nwk::NwkStyle;
use treetime_utils::io::file::is_path_stdout;
use treetime_utils::make_error;

pub fn plan(request: &OutputPlanRequest) -> Result<ResolvedOutputs, Report> {
  let command = request.command;
  let stem = command.stem();
  let partitioned = request.partition_count.get() > 1;
  let styles = effective_nwk_styles(command, &request.nwk_styles);

  let mut tree_outputs: BTreeMap<TreeWriteKind, PathBuf> = BTreeMap::new();
  let mut non_tree_outputs: BTreeMap<OutputSelection, PlannedFile> = BTreeMap::new();

  let mut overridden_tree: BTreeSet<OutputSelection> = BTreeSet::new();
  for (&variant, path) in &request.tree_overrides {
    overridden_tree.insert(variant);
    if variant.is_styled_tree() {
      for (style, p) in expand_override_styles(path, &styles) {
        tree_outputs.insert(styled_tree_write_kind(variant, style), p);
      }
    } else if let Some(kind) = variant.to_tree_write_kind() {
      tree_outputs.insert(kind, path.clone());
    }
  }

  for (&sel, path) in &request.non_tree_overrides {
    non_tree_outputs.insert(
      sel,
      PlannedFile {
        path: path.clone(),
        requested: Requested::Named,
      },
    );
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
          for (style, p) in expand_outputall_styles(dir, stem, partitioned, command, variant, &styles) {
            tree_outputs.entry(styled_tree_write_kind(variant, style)).or_insert(p);
          }
        } else if let Some(kind) = variant.to_tree_write_kind() {
          let stem = default_stem(stem, partitioned && command.is_per_partition_tree(kind));
          tree_outputs
            .entry(kind)
            .or_insert_with(|| dir.join(format!("{stem}{}", variant.extension())));
        }
      } else {
        let stem = default_stem(stem, partitioned && command.is_per_partition(variant));
        non_tree_outputs.entry(variant).or_insert_with(|| PlannedFile {
          path: dir.join(format!("{stem}{}", variant.extension())),
          requested: Requested::All,
        });
      }
    }
  } else if !request.selection.is_empty() {
    return Err(Report::new(MissingOutputs::SelectionWithoutFolder));
  }

  if tree_outputs.is_empty() && non_tree_outputs.is_empty() {
    return Err(Report::new(MissingOutputs::NoOutputs));
  }

  let resolved = expand_partitions(command, request.partition_count, tree_outputs, non_tree_outputs)?;
  ensure_unique_destinations(
    resolved
      .planned_paths()
      .into_iter()
      .map(|planned| (planned.label, planned.path)),
  )?;
  Ok(resolved)
}

fn expand_partitions(
  command: CommandKind,
  partition_count: NonZeroUsize,
  tree_outputs: BTreeMap<TreeWriteKind, PathBuf>,
  non_tree_outputs: BTreeMap<OutputSelection, PlannedFile>,
) -> Result<ResolvedOutputs, Report> {
  let partitions = (1..=partition_count.get())
    .map(|partition| partition.to_string())
    .collect_vec();
  let mut resolved = ResolvedOutputs {
    partition_count,
    tree_outputs: BTreeMap::new(),
    non_tree_outputs: BTreeMap::new(),
    tree_outputs_by_partition: BTreeMap::new(),
    non_tree_outputs_by_partition: BTreeMap::new(),
  };
  for (kind, path) in tree_outputs {
    let label = tree_write_kind_label(kind);
    if command.is_per_partition_tree(kind)
      || (command.reads_alignments() && PathPlaceholder::Partition.is_in_file_name(&path))
    {
      let paths = expand_path_template(&label, &path, PathPlaceholder::Partition, &partitions)?;
      resolved.tree_outputs_by_partition.insert(kind, paths);
    } else {
      ensure_no_partition_placeholder(&label, &path)?;
      resolved.tree_outputs.insert(kind, path);
    }
  }
  for (selection, file) in non_tree_outputs {
    let label = selection.flag_name();
    if command.is_per_partition(selection) {
      let files = expand_path_template(label, &file.path, PathPlaceholder::Partition, &partitions)?
        .into_iter()
        .map(|path| PlannedFile {
          path,
          requested: file.requested,
        })
        .collect();
      resolved.non_tree_outputs_by_partition.insert(selection, files);
    } else {
      ensure_no_partition_placeholder(label, &file.path)?;
      resolved.non_tree_outputs.insert(selection, file);
    }
  }
  Ok(resolved)
}

pub fn expand_path_template(
  flag: &str,
  template: &Path,
  placeholder: PathPlaceholder,
  values: &[String],
) -> Result<Vec<PathBuf>, Report> {
  let Some(name) = placeholder.file_name_template(template) else {
    if values.len() > 1 {
      return make_error!(
        "{flag} needs a {noun} placeholder ({tokens}) in its file name when {reason}, \
         otherwise each {noun} overwrites the same output file",
        noun = placeholder.noun(),
        tokens = placeholder.quoted_tokens(),
        reason = placeholder.several_values_reason(),
      );
    }
    return Ok(values.iter().map(|_| template.to_path_buf()).collect());
  };
  values
    .iter()
    .map(|value| {
      if value.chars().any(is_separator) {
        return make_error!(
          "{} '{value}' contains a path separator and cannot be part of the file name of {flag}",
          placeholder.noun()
        );
      }
      Ok(template.with_file_name(placeholder.substitute(name, value)))
    })
    .collect()
}

fn ensure_no_partition_placeholder(label: &str, path: &Path) -> Result<(), Report> {
  if PathPlaceholder::Partition.is_in_file_name(path) {
    return make_error!(
      "{label} is written once per run, so its file name cannot contain {}",
      PathPlaceholder::Partition.quoted_tokens()
    );
  }
  Ok(())
}

fn tree_write_kind_label(kind: TreeWriteKind) -> Cow<'static, str> {
  let flag = tree_write_kind_selection(kind).flag_name();
  match kind {
    TreeWriteKind::Nwk(style) | TreeWriteKind::Nexus(style) => {
      Cow::Owned(format!("{flag} (style {})", nwk_style_name(style)))
    },
    TreeWriteKind::Auspice
    | TreeWriteKind::MatPb
    | TreeWriteKind::MatJson
    | TreeWriteKind::GraphJson
    | TreeWriteKind::Dot => Cow::Borrowed(flag),
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
  partitioned: bool,
  command: CommandKind,
  variant: OutputSelection,
  styles: &[NwkStyle],
) -> Vec<(NwkStyle, PathBuf)> {
  let multi = styles.len() > 1;
  let primary = variant.extension();
  styles
    .iter()
    .map(|&style| {
      let secondary = if multi { nwk_style_secondary_ext(style) } else { "" };
      let stem = default_stem(
        stem,
        partitioned && command.is_per_partition_tree(styled_tree_write_kind(variant, style)),
      );
      (style, dir.join(format!("{stem}{secondary}{primary}")))
    })
    .collect()
}

fn default_stem(stem: &str, per_partition: bool) -> Cow<'_, str> {
  match (per_partition, PathPlaceholder::Partition.tokens()) {
    (true, [token, ..]) => Cow::Owned(format!("{stem}.{token}")),
    _ => Cow::Borrowed(stem),
  }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PathPlaceholder {
  Partition,
  Cds,
}

impl PathPlaceholder {
  pub const fn tokens(self) -> &'static [&'static str] {
    match self {
      Self::Partition => &["{partition}"],
      Self::Cds => &["{cds}", "%GENE"],
    }
  }

  pub fn is_in(self, text: &str) -> bool {
    self.tokens().iter().any(|token| text.contains(token))
  }

  pub fn is_in_file_name(self, path: &Path) -> bool {
    self.file_name_template(path).is_some()
  }

  pub fn substitute(self, text: &str, value: &str) -> String {
    self
      .tokens()
      .iter()
      .fold(text.to_owned(), |expanded, token| expanded.replace(token, value))
  }

  pub fn quoted_tokens(self) -> String {
    self.tokens().iter().map(|token| format!("'{token}'")).join(" or ")
  }

  fn file_name_template(self, path: &Path) -> Option<&str> {
    path
      .file_name()
      .and_then(|name| name.to_str())
      .filter(|name| self.is_in(name))
  }

  const fn noun(self) -> &'static str {
    match self {
      Self::Partition => "partition",
      Self::Cds => "CDS",
    }
  }

  const fn several_values_reason(self) -> &'static str {
    match self {
      Self::Partition => "the run has several alignment files",
      Self::Cds => "reconstructing several CDSes",
    }
  }
}

fn styled_tree_write_kind(variant: OutputSelection, style: NwkStyle) -> TreeWriteKind {
  match variant {
    OutputSelection::Nwk => TreeWriteKind::Nwk(style),
    OutputSelection::Nexus => TreeWriteKind::Nexus(style),
    _ => unreachable!("styled_tree_write_kind called on non-styled variant"),
  }
}

fn nwk_style_name(style: NwkStyle) -> &'static str {
  match style {
    NwkStyle::Plain => "plain",
    NwkStyle::Beast => "beast",
    NwkStyle::Nhx => "nhx",
  }
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

  pub partition_count: NonZeroUsize,

  pub output_all: Option<PathBuf>,

  pub nwk_styles: Vec<NwkStyle>,

  pub selection: Vec<OutputSelection>,

  pub tree_overrides: BTreeMap<OutputSelection, PathBuf>,

  pub non_tree_overrides: BTreeMap<OutputSelection, PathBuf>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Display, Error)]
pub enum MissingOutputs {
  #[display("--output-selection requires --output-all")]
  SelectionWithoutFolder,
  #[display(
    "No output flags provided. At least one is required: --output-all or one of the --output-tree-* / --output-* flags"
  )]
  NoOutputs,
}

pub struct ResolvedOutputs {
  pub partition_count: NonZeroUsize,
  pub tree_outputs: BTreeMap<TreeWriteKind, PathBuf>,
  pub non_tree_outputs: BTreeMap<OutputSelection, PlannedFile>,
  pub tree_outputs_by_partition: BTreeMap<TreeWriteKind, Vec<PathBuf>>,
  pub non_tree_outputs_by_partition: BTreeMap<OutputSelection, Vec<PlannedFile>>,
}

impl ResolvedOutputs {
  pub fn paths_by_selection(&self) -> BTreeMap<OutputSelection, Vec<PathBuf>> {
    let mut by_selection: BTreeMap<OutputSelection, Vec<PathBuf>> = BTreeMap::new();
    for planned in self.planned_paths() {
      by_selection
        .entry(planned.selection)
        .or_default()
        .push(planned.path.to_path_buf());
    }
    for paths in by_selection.values_mut() {
      paths.sort();
    }
    by_selection
  }

  pub fn partitions(&self) -> RangeInclusive<usize> {
    1..=self.partition_count.get()
  }

  pub fn path(&self, selection: OutputSelection) -> Option<&Path> {
    self.non_tree_outputs.get(&selection).map(|file| file.path.as_path())
  }

  pub fn partition_files(&self, selection: OutputSelection) -> &[PlannedFile] {
    self
      .non_tree_outputs_by_partition
      .get(&selection)
      .map_or(&[], Vec::as_slice)
  }

  pub fn partition_paths(&self, selection: OutputSelection) -> BTreeMap<usize, PathBuf> {
    self
      .partition_files(selection)
      .iter()
      .zip(1..)
      .map(|(file, partition)| (partition, file.path.clone()))
      .collect()
  }

  pub fn partition_path(&self, selection: OutputSelection, partition: usize) -> Option<&Path> {
    partition_entry(self.partition_files(selection), partition).map(|file| file.path.as_path())
  }

  pub fn tree_partition_path(&self, kind: TreeWriteKind, partition: usize) -> Option<&Path> {
    let paths = self.tree_outputs_by_partition.get(&kind).map_or(&[][..], Vec::as_slice);
    partition_entry(paths, partition).map(PathBuf::as_path)
  }

  pub fn partition_tree_outputs(&self, partition: usize) -> BTreeMap<TreeWriteKind, PathBuf> {
    self
      .tree_outputs_by_partition
      .keys()
      .filter_map(|&kind| Some((kind, self.tree_partition_path(kind, partition)?.to_path_buf())))
      .collect()
  }

  pub fn tree_write_kinds(&self) -> impl Iterator<Item = TreeWriteKind> + '_ {
    self
      .tree_outputs
      .keys()
      .chain(self.tree_outputs_by_partition.keys())
      .copied()
  }

  pub fn tree_based_paths(&self) -> Vec<&Path> {
    self
      .planned_paths()
      .into_iter()
      .filter(|planned| planned.selection.reads_tree())
      .map(|planned| planned.path)
      .collect()
  }

  pub fn ensure_unique_with_expansion<'a>(
    &'a self,
    selection: OutputSelection,
    expanded: impl IntoIterator<Item = (String, &'a Path)>,
  ) -> Result<(), Report> {
    let planned = self
      .planned_paths()
      .into_iter()
      .filter(|planned| planned.selection != selection)
      .map(|planned| (planned.label, planned.path));
    let expanded = expanded.into_iter().map(|(label, path)| (Cow::Owned(label), path));
    ensure_unique_destinations(planned.chain(expanded))
  }

  fn planned_paths(&self) -> Vec<PlannedPath<'_>> {
    let run_trees = self
      .tree_outputs
      .iter()
      .map(|(kind, path)| PlannedPath::run(tree_write_kind_selection(*kind), path));
    let run_files = self
      .non_tree_outputs
      .iter()
      .map(|(selection, file)| PlannedPath::run(*selection, &file.path));
    let partition_trees = self.tree_outputs_by_partition.iter().flat_map(|(kind, paths)| {
      PlannedPath::partitions(tree_write_kind_selection(*kind), paths.iter().map(PathBuf::as_path))
    });
    let partition_files = self
      .non_tree_outputs_by_partition
      .iter()
      .flat_map(|(selection, files)| PlannedPath::partitions(*selection, files.iter().map(|file| file.path.as_path())));
    run_trees
      .chain(run_files)
      .chain(partition_trees)
      .chain(partition_files)
      .collect()
  }
}

struct PlannedPath<'a> {
  selection: OutputSelection,
  label: Cow<'static, str>,
  path: &'a Path,
}

impl<'a> PlannedPath<'a> {
  fn run(selection: OutputSelection, path: &'a Path) -> Self {
    Self {
      selection,
      label: Cow::Borrowed(selection.flag_name()),
      path,
    }
  }

  fn partitions(selection: OutputSelection, paths: impl ExactSizeIterator<Item = &'a Path>) -> Vec<Self> {
    if paths.len() == 1 {
      return paths.map(|path| Self::run(selection, path)).collect();
    }
    paths
      .zip(1..)
      .map(|(path, partition)| Self {
        selection,
        label: Cow::Owned(format!("{} for partition {partition}", selection.flag_name())),
        path,
      })
      .collect()
  }
}

fn partition_entry<T>(entries: &[T], partition: usize) -> Option<&T> {
  entries.get(partition.checked_sub(1)?)
}

pub fn output_unavailable(
  selection: OutputSelection,
  file: &PlannedFile,
  reason: &str,
  log: &dyn LogSink,
) -> Result<(), Report> {
  match file.requested {
    Requested::Named => make_error!("{} was requested, but {reason}", selection.flag_name()),
    Requested::All => {
      progress_debug!(log, "Not writing '{}': {reason}", file.path.display());
      Ok(())
    },
  }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PlannedFile {
  pub path: PathBuf,
  pub requested: Requested,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Requested {
  Named,
  All,
}

fn tree_write_kind_selection(kind: TreeWriteKind) -> OutputSelection {
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
  Homoplasy,
}

impl CommandKind {
  pub fn all_selectable(self) -> BTreeSet<OutputSelection> {
    &Self::available_tree_outputs() | &self.non_tree_outputs()
  }

  pub fn default_outputs(self) -> BTreeSet<OutputSelection> {
    let mut non_tree = self.non_tree_outputs();
    let tree_defaults = if self == Self::Timetree {
      non_tree.remove(&OutputSelection::ClockCsv);
      btreeset![OutputSelection::Nwk, OutputSelection::Nexus, OutputSelection::Auspice]
    } else {
      btreeset![OutputSelection::Nwk, OutputSelection::Nexus]
    };
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
        OutputSelection::ClockCsv,
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
      Self::Clock => btreeset![
        OutputSelection::ClockModel,
        OutputSelection::ClockCsv,
        OutputSelection::ClockChartSvg,
        OutputSelection::ClockChartPng,
        OutputSelection::ClockChartPdf
      ],
      Self::Prune => btreeset![OutputSelection::Gtr],
      Self::Homoplasy => btreeset![OutputSelection::HomoplasyStats, OutputSelection::HomoplasyReport],
    }
  }

  pub(crate) fn stem(self) -> &'static str {
    match self {
      Self::Ancestral => "ancestral",
      Self::Timetree => "timetree",
      Self::Optimize => "optimize",
      Self::Mugration => "mugration",
      Self::Clock => "clock",
      Self::Prune => "prune",
      Self::Homoplasy => "homoplasy",
    }
  }

  fn reads_alignments(self) -> bool {
    match self {
      Self::Ancestral | Self::Timetree | Self::Optimize | Self::Prune | Self::Homoplasy => true,
      Self::Mugration | Self::Clock => false,
    }
  }

  fn is_per_partition(self, selection: OutputSelection) -> bool {
    match selection {
      OutputSelection::ReconstructedNucFasta => true,
      OutputSelection::AugurNodeData => self == Self::Ancestral,
      _ => false,
    }
  }

  fn is_per_partition_tree(self, kind: TreeWriteKind) -> bool {
    self.reads_alignments() && kind.is_per_partition()
  }

  fn default_nwk_styles(self) -> Vec<NwkStyle> {
    vec![NwkStyle::Plain]
  }
}

fn ensure_unique_destinations<'a>(
  outputs: impl IntoIterator<Item = (Cow<'static, str>, &'a Path)>,
) -> Result<(), Report> {
  let mut destinations: BTreeMap<&Path, Cow<'static, str>> = BTreeMap::new();
  for (label, path) in outputs {
    let destination = if is_path_stdout(path) { Path::new("-") } else { path };
    if let Some(previous) = destinations.get(destination) {
      return make_error!(
        "Output destination '{}' is selected more than once ({previous} and {label})",
        path.display()
      );
    }
    destinations.insert(destination, label);
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
  JsonSchema,
  AsRefStr,
  EnumString,
  EnumIter,
  Serialize,
  Deserialize,
)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
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
  ClockChartSvg,
  ClockChartPng,
  ClockChartPdf,
  HomoplasyStats,
  HomoplasyReport,
}

impl OutputSelection {
  pub const fn description(self) -> &'static str {
    self.facts().description
  }

  pub const fn table_format(self) -> Option<TableFormat> {
    self.facts().table_format
  }

  pub const fn flag_name(self) -> &'static str {
    self.facts().flag_name
  }

  fn is_tree(self) -> bool {
    matches!(
      self,
      Self::Nwk | Self::Nexus | Self::Auspice | Self::MatPb | Self::MatJson | Self::GraphJson | Self::Dot
    )
  }

  const fn reads_tree(self) -> bool {
    self.facts().reads_tree
  }

  fn is_styled_tree(self) -> bool {
    matches!(self, Self::Nwk | Self::Nexus)
  }

  fn is_meta(self) -> bool {
    matches!(self, Self::All)
  }

  const fn extension(self) -> &'static str {
    self.facts().extension
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

  #[rustfmt::skip]
  const fn facts(self) -> OutputFacts {
    const CSV: Option<TableFormat> = Some(TableFormat::Csv);
    const TSV: Option<TableFormat> = Some(TableFormat::Tsv);
    match self {
      Self::All                   => OutputFacts::new("",                              "--output-selection=all",          None, false, "Every output of the command"),
      Self::Nwk                   => OutputFacts::new(".nwk",                          "--output-tree-nwk",               None, true,  "Tree in Newick format"),
      Self::Nexus                 => OutputFacts::new(".nexus",                        "--output-tree-nexus",             None, true,  "Tree with node annotations in Nexus format"),
      Self::Auspice               => OutputFacts::new(".auspice.json",                 "--output-tree-auspice",           None, true,  "Tree for Auspice and Nextstrain"),
      Self::MatPb                 => OutputFacts::new(".mat.pb",                       "--output-tree-mat-pb",            None, true,  "Mutation-annotated tree (UShER protobuf)"),
      Self::MatJson               => OutputFacts::new(".mat.json",                     "--output-tree-mat-json",          None, true,  "Mutation-annotated tree (JSON)"),
      Self::GraphJson             => OutputFacts::new(".graph.json",                   "--output-tree-graph-json",        None, false, "Tree as a graph (JSON)"),
      Self::Dot                   => OutputFacts::new(".dot",                          "--output-tree-dot",               None, false, "Tree as a graph (Graphviz)"),
      Self::AugurNodeData         => OutputFacts::new(".augur-node-data.json",         "--output-augur-node-data",        None, true,  "Node data for augur: dates, intervals, mutations"),
      Self::Gtr                   => OutputFacts::new(".gtr.json",                     "--output-gtr",                    None, false, "Substitution model"),
      Self::ClockModel            => OutputFacts::new(".clock-model.json",             "--output-clock-model",            None, false, "Clock rate, intercept and regression statistics"),
      Self::ConfidenceTsv         => OutputFacts::new(".confidence.tsv",               "--output-confidence-tsv",         TSV,  false, "Date intervals of every node"),
      Self::ConfidenceCsv         => OutputFacts::new(".confidence.csv",               "--output-confidence-csv",         CSV,  true,  "State probabilities of every node"),
      Self::ReconstructedNucFasta => OutputFacts::new(".reconstructed-nuc.fasta",      "--output-reconstructed-nuc-fasta", None, false, "Sequences of samples and ancestors"),
      Self::ReconstructedAaFasta  => OutputFacts::new(".reconstructed-aa.{cds}.fasta", "--output-reconstructed-aa-fasta", None, false, "Protein sequences of samples and ancestors"),
      Self::TraitsCsv             => OutputFacts::new(".traits.csv",                   "--output-traits-csv",             CSV,  true,  "Inferred state of every node"),
      Self::ClockCsv              => OutputFacts::new(".clock.csv",                    "--output-clock-csv",              CSV,  false, "Root-to-tip distance, date and clock prediction of every node"),
      Self::Tracelog              => OutputFacts::new(".tracelog.csv",                 "--output-tracelog",               CSV,  false, "Convergence values of every iteration"),
      Self::CoalescentTsv         => OutputFacts::new(".coalescent.tsv",               "--output-coalescent-tsv",         TSV,  false, "Coalescent time scale and effective population size"),
      Self::CoalescentCsv         => OutputFacts::new(".coalescent.csv",               "--output-coalescent-csv",         CSV,  false, "Coalescent time scale and effective population size"),
      Self::CoalescentJson        => OutputFacts::new(".coalescent.json",              "--output-coalescent-json",        None, false, "Coalescent model and its likelihood"),
      Self::ClockChartSvg         => OutputFacts::new(".svg",                          "--output-clock-chart-svg",        None, false, "Root-to-tip regression chart (SVG)"),
      Self::ClockChartPng         => OutputFacts::new(".png",                          "--output-clock-chart-png",        None, false, "Root-to-tip regression chart (PNG)"),
      Self::ClockChartPdf         => OutputFacts::new(".pdf",                          "--output-clock-chart-pdf",        None, false, "Root-to-tip regression chart (PDF)"),
      Self::HomoplasyStats        => OutputFacts::new(".stats.json",                   "--output-homoplasy-stats",        None, false, "Recurrent mutations: counts, histograms and ranked lists (JSON)"),
      Self::HomoplasyReport       => OutputFacts::new(".report.txt",                   "--output-homoplasy-report",       None, false, "Recurrent mutations: text report"),
    }
  }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum TreeWriteKind {
  Nwk(NwkStyle),
  Nexus(NwkStyle),
  Auspice,
  MatPb,
  MatJson,
  GraphJson,
  Dot,
}

impl TreeWriteKind {
  pub fn is_per_partition(self) -> bool {
    match self {
      Self::Nwk(style) | Self::Nexus(style) => style != NwkStyle::Plain,
      Self::Auspice | Self::MatPb | Self::MatJson => true,
      Self::GraphJson | Self::Dot => false,
    }
  }
}

struct OutputFacts {
  extension: &'static str,
  flag_name: &'static str,
  table_format: Option<TableFormat>,
  reads_tree: bool,
  description: &'static str,
}

impl OutputFacts {
  const fn new(
    extension: &'static str,
    flag_name: &'static str,
    table_format: Option<TableFormat>,
    reads_tree: bool,
    description: &'static str,
  ) -> Self {
    Self {
      extension,
      flag_name,
      table_format,
      reads_tree,
      description,
    }
  }
}
