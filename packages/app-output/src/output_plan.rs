use crate::command_outputs::CommandKind;
use crate::output_conditions::{OutputCondition, OutputSettings, UnavailableOutput, UnmetCondition};
use crate::output_format::FileFormat;
use derive_more::{Display, Error};
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};
use std::num::NonZeroUsize;
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf, is_separator};
use strum_macros::{AsRefStr, EnumIter, EnumString};
use treetime::progress::LogSink;
use treetime::{progress_debug, progress_info};
use treetime_io::nwk::NwkStyle;
use treetime_utils::io::file::is_path_stdout;
use treetime_utils::make_error;

pub fn plan(request: &OutputPlanRequest) -> Result<ResolvedOutputs, Report> {
  let command = request.command;
  let settings = request.settings;
  let stem = command.stem();
  let per_alignment = request.alignment_count.get() > 1;

  let mut tree_outputs: BTreeMap<TreeWriteKind, PathBuf> = BTreeMap::new();
  let mut non_tree_outputs: BTreeMap<OutputSelection, PlannedFile> = BTreeMap::new();

  let mut overridden_tree: BTreeSet<OutputSelection> = BTreeSet::new();
  for (&variant, path) in &request.tree_overrides {
    overridden_tree.insert(variant);
    if let Some(kind) = variant.to_tree_write_kind(request.newick_annotations) {
      tree_outputs.insert(kind, path.clone());
    }
  }

  for (&sel, path) in &request.non_tree_overrides {
    ensure_available(settings, sel, Requested::Named)?;
    non_tree_outputs.insert(
      sel,
      PlannedFile::new(path.clone(), Requested::Named, settings.condition(sel)),
    );
  }

  if let Some(dir) = &request.output_all {
    let mut effective: BTreeSet<OutputSelection> = if request.selection.is_empty() {
      command.default_outputs(settings)
    } else if request.selection.contains(&OutputSelection::All) {
      command.all_selectable(settings)
    } else {
      for &selection in &request.selection {
        ensure_available(settings, selection, Requested::All)?;
      }
      request.selection.iter().copied().collect()
    };
    effective.extend(command.added_outputs(settings).map(|(selection, _)| selection));

    for &variant in &effective {
      if variant.is_meta() {
        continue;
      }
      if variant.is_tree() {
        if overridden_tree.contains(&variant) {
          continue;
        }
        if let Some(kind) = variant.to_tree_write_kind(request.newick_annotations) {
          let stem = default_stem(stem, per_alignment && command.is_per_alignment_tree(kind));
          tree_outputs
            .entry(kind)
            .or_insert_with(|| dir.join(format!("{stem}{}", variant.extension())));
        }
      } else {
        let stem = default_stem(stem, per_alignment && command.is_per_alignment(variant));
        non_tree_outputs.entry(variant).or_insert_with(|| {
          PlannedFile::new(
            dir.join(format!("{stem}{}", variant.extension())),
            Requested::All,
            settings.condition(variant),
          )
        });
      }
    }
  } else if !request.selection.is_empty() {
    return Err(Report::new(MissingOutputs::SelectionWithoutFolder));
  }

  if request.output_all.is_none() {
    for (selection, setting) in command.added_outputs(settings) {
      if !non_tree_outputs.contains_key(&selection) {
        return Err(Report::new(AddedOutputWithoutFolder { selection, setting }));
      }
    }
  }
  if tree_outputs.is_empty() && non_tree_outputs.is_empty() {
    return Err(Report::new(MissingOutputs::NoOutputs));
  }

  let resolved = expand_alignments(command, request.alignment_count, tree_outputs, non_tree_outputs)?;
  for planned in resolved.planned_paths() {
    planned.selection.file_format().check(planned.selection, planned.path)?;
  }
  ensure_unique_destinations(
    resolved
      .planned_paths()
      .into_iter()
      .map(|planned| (planned.label, planned.path)),
  )?;
  Ok(resolved)
}

fn expand_alignments(
  command: CommandKind,
  alignment_count: NonZeroUsize,
  tree_outputs: BTreeMap<TreeWriteKind, PathBuf>,
  non_tree_outputs: BTreeMap<OutputSelection, PlannedFile>,
) -> Result<ResolvedOutputs, Report> {
  let alignments = (1..=alignment_count.get())
    .map(|alignment| alignment.to_string())
    .collect_vec();
  let mut resolved = ResolvedOutputs {
    alignment_count,
    tree_outputs: BTreeMap::new(),
    non_tree_outputs: BTreeMap::new(),
    tree_outputs_by_alignment: BTreeMap::new(),
    non_tree_outputs_by_alignment: BTreeMap::new(),
  };
  for (kind, path) in tree_outputs {
    let label = tree_write_kind_selection(kind).flag_name();
    if command.is_per_alignment_tree(kind)
      || (command.reads_alignments() && PathPlaceholder::Alignment.is_in_file_name(&path))
    {
      let paths = expand_path_template(label, &path, PathPlaceholder::Alignment, &alignments)?;
      resolved.tree_outputs_by_alignment.insert(kind, paths);
    } else {
      ensure_no_alignment_placeholder(label, &path)?;
      resolved.tree_outputs.insert(kind, path);
    }
  }
  for (selection, file) in non_tree_outputs {
    let label = selection.flag_name();
    if command.is_per_alignment(selection) {
      let files = expand_path_template(label, &file.path, PathPlaceholder::Alignment, &alignments)?
        .into_iter()
        .map(|path| PlannedFile::new(path, file.requested, file.condition))
        .collect();
      resolved.non_tree_outputs_by_alignment.insert(selection, files);
    } else {
      ensure_no_alignment_placeholder(label, &file.path)?;
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
        "{flag} needs the placeholder {tokens} in its file name when {reason}, \
         otherwise the output of each {noun} overwrites the same file",
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

fn ensure_no_alignment_placeholder(label: &str, path: &Path) -> Result<(), Report> {
  if PathPlaceholder::Alignment.is_in_file_name(path) {
    return make_error!(
      "{label} is written once per run, so its file name cannot contain {}",
      PathPlaceholder::Alignment.quoted_tokens()
    );
  }
  Ok(())
}

fn default_stem(stem: &str, per_alignment: bool) -> Cow<'_, str> {
  match (per_alignment, PathPlaceholder::Alignment.tokens()) {
    (true, [token, ..]) => Cow::Owned(format!("{stem}.{token}")),
    _ => Cow::Borrowed(stem),
  }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PathPlaceholder {
  Alignment,
  Cds,
}

impl PathPlaceholder {
  pub const fn tokens(self) -> &'static [&'static str] {
    match self {
      Self::Alignment => &["{alignment}"],
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
      Self::Alignment => "alignment file",
      Self::Cds => "CDS",
    }
  }

  const fn several_values_reason(self) -> &'static str {
    match self {
      Self::Alignment => "the run has several alignment files",
      Self::Cds => "reconstructing several CDSes",
    }
  }
}

pub struct OutputPlanRequest {
  pub command: CommandKind,

  pub alignment_count: NonZeroUsize,

  pub output_all: Option<PathBuf>,

  pub newick_annotations: NwkStyle,

  pub selection: Vec<OutputSelection>,

  pub tree_overrides: BTreeMap<OutputSelection, PathBuf>,

  pub non_tree_overrides: BTreeMap<OutputSelection, PathBuf>,

  pub settings: OutputSettings,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Display, Error)]
pub enum MissingOutputs {
  #[display("--output-selection requires --output-all")]
  SelectionWithoutFolder,
  #[display(
    "No outputs selected: add -O <DIR> to write the default outputs into a folder, or a per-file --output-* flag"
  )]
  NoOutputs,
}

pub struct ResolvedOutputs {
  pub alignment_count: NonZeroUsize,
  pub tree_outputs: BTreeMap<TreeWriteKind, PathBuf>,
  pub non_tree_outputs: BTreeMap<OutputSelection, PlannedFile>,
  pub tree_outputs_by_alignment: BTreeMap<TreeWriteKind, Vec<PathBuf>>,
  pub non_tree_outputs_by_alignment: BTreeMap<OutputSelection, Vec<PlannedFile>>,
}

impl ResolvedOutputs {
  #[expect(
    clippy::zero_sized_map_values,
    reason = "the map names the condition of each output, of which there is one kind"
  )]
  pub fn conditional_selections(&self) -> BTreeMap<OutputSelection, OutputCondition> {
    self
      .non_tree_outputs
      .iter()
      .map(|(&selection, file)| (selection, file.condition))
      .chain(
        self
          .non_tree_outputs_by_alignment
          .iter()
          .flat_map(|(&selection, files)| files.iter().map(move |file| (selection, file.condition))),
      )
      .filter_map(|(selection, condition)| Some((selection, condition?)))
      .collect()
  }

  pub fn settle(&mut self, condition: OutputCondition, holds: bool, log: &dyn LogSink) -> Result<(), UnmetCondition> {
    let files = self
      .non_tree_outputs
      .iter_mut()
      .map(|(&selection, file)| (selection, file))
      .chain(
        self
          .non_tree_outputs_by_alignment
          .iter_mut()
          .flat_map(|(&selection, files)| files.iter_mut().map(move |file| (selection, file))),
      )
      .filter(|(_, file)| file.condition == Some(condition));
    let mut skipped = Vec::new();
    for (selection, file) in files {
      if holds {
        file.condition = None;
      } else if file.requested == Requested::Named {
        return Err(UnmetCondition { selection, condition });
      } else {
        progress_info!(
          log,
          "Not writing '{}': {}",
          file.path.display(),
          condition.unmet_reason()
        );
        skipped.push(selection);
      }
    }
    for selection in skipped {
      self.non_tree_outputs.remove(&selection);
      self.non_tree_outputs_by_alignment.remove(&selection);
    }
    Ok(())
  }

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

  pub fn alignments(&self) -> RangeInclusive<usize> {
    1..=self.alignment_count.get()
  }

  pub fn path(&self, selection: OutputSelection) -> Option<&Path> {
    self.non_tree_outputs.get(&selection).map(|file| file.path.as_path())
  }

  pub fn alignment_files(&self, selection: OutputSelection) -> &[PlannedFile] {
    self
      .non_tree_outputs_by_alignment
      .get(&selection)
      .map_or(&[], Vec::as_slice)
  }

  pub fn alignment_paths(&self, selection: OutputSelection) -> BTreeMap<usize, PathBuf> {
    self
      .alignment_files(selection)
      .iter()
      .zip(1..)
      .map(|(file, alignment)| (alignment, file.path.clone()))
      .collect()
  }

  pub fn alignment_path(&self, selection: OutputSelection, alignment: usize) -> Option<&Path> {
    alignment_entry(self.alignment_files(selection), alignment).map(|file| file.path.as_path())
  }

  pub fn tree_alignment_path(&self, kind: TreeWriteKind, alignment: usize) -> Option<&Path> {
    let paths = self.tree_outputs_by_alignment.get(&kind).map_or(&[][..], Vec::as_slice);
    alignment_entry(paths, alignment).map(PathBuf::as_path)
  }

  pub fn alignment_tree_outputs(&self, alignment: usize) -> BTreeMap<TreeWriteKind, PathBuf> {
    self
      .tree_outputs_by_alignment
      .keys()
      .filter_map(|&kind| Some((kind, self.tree_alignment_path(kind, alignment)?.to_path_buf())))
      .collect()
  }

  pub fn tree_write_kinds(&self) -> impl Iterator<Item = TreeWriteKind> + '_ {
    self
      .tree_outputs
      .keys()
      .chain(self.tree_outputs_by_alignment.keys())
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
    let alignment_trees = self.tree_outputs_by_alignment.iter().flat_map(|(kind, paths)| {
      PlannedPath::per_alignment(tree_write_kind_selection(*kind), paths.iter().map(PathBuf::as_path))
    });
    let alignment_files = self
      .non_tree_outputs_by_alignment
      .iter()
      .flat_map(|(selection, files)| {
        PlannedPath::per_alignment(*selection, files.iter().map(|file| file.path.as_path()))
      });
    run_trees
      .chain(run_files)
      .chain(alignment_trees)
      .chain(alignment_files)
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

  fn per_alignment(selection: OutputSelection, paths: impl ExactSizeIterator<Item = &'a Path>) -> Vec<Self> {
    if paths.len() == 1 {
      return paths.map(|path| Self::run(selection, path)).collect();
    }
    paths
      .zip(1..)
      .map(|(path, alignment)| Self {
        selection,
        label: Cow::Owned(format!("{} for alignment {alignment}", selection.flag_name())),
        path,
      })
      .collect()
  }
}

fn alignment_entry<T>(entries: &[T], alignment: usize) -> Option<&T> {
  entries.get(alignment.checked_sub(1)?)
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
  pub condition: Option<OutputCondition>,
}

impl PlannedFile {
  pub const fn new(path: PathBuf, requested: Requested, condition: Option<OutputCondition>) -> Self {
    Self {
      path,
      requested,
      condition,
    }
  }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Requested {
  Named,
  All,
}

fn ensure_available(
  settings: OutputSettings,
  selection: OutputSelection,
  requested: Requested,
) -> Result<(), UnavailableOutput> {
  match settings.unavailable(selection) {
    Some(cause) => Err(UnavailableOutput {
      selection,
      requested,
      cause,
    }),
    None => Ok(()),
  }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Display, Error)]
#[display(
  "--{} has no effect without --output-all: add --output-all <DIR>, or give the file with {} <FILE>",
  setting.replace('_', "-"),
  selection.flag_name()
)]
pub struct AddedOutputWithoutFolder {
  pub selection: OutputSelection,
  pub setting: &'static str,
}

fn tree_write_kind_selection(kind: TreeWriteKind) -> OutputSelection {
  match kind {
    TreeWriteKind::Nwk(_) => OutputSelection::Nwk,
    TreeWriteKind::Nexus => OutputSelection::Nexus,
    TreeWriteKind::Auspice => OutputSelection::Auspice,
    TreeWriteKind::MatPb => OutputSelection::MatPb,
    TreeWriteKind::MatJson => OutputSelection::MatJson,
    TreeWriteKind::GraphJson => OutputSelection::GraphJson,
    TreeWriteKind::Dot => OutputSelection::Dot,
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
/// Tree variants do not encode the separately selected Newick annotations.
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

  NodeData,
  Model,
  ClockModel,
  Dates,
  TraitProbabilities,
  Sequences,
  Translations,
  Traits,
  Rtt,
  Tracelog,
  Coalescent,
  CoalescentJson,
  RttChart,
  HomoplasyStats,
  HomoplasyReport,
}

impl OutputSelection {
  pub const fn description(self) -> &'static str {
    self.facts().description
  }

  pub const fn file_format(self) -> FileFormat {
    self.facts().format
  }

  pub const fn flag_name(self) -> &'static str {
    self.facts().flag_name
  }

  pub fn setting_key(self) -> String {
    self.flag_name().trim_start_matches("--").replace('-', "_")
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

  fn is_meta(self) -> bool {
    matches!(self, Self::All)
  }

  const fn extension(self) -> &'static str {
    self.facts().extension
  }

  fn to_tree_write_kind(self, newick_annotations: NwkStyle) -> Option<TreeWriteKind> {
    match self {
      Self::Nwk => Some(TreeWriteKind::Nwk(newick_annotations)),
      Self::Nexus => Some(TreeWriteKind::Nexus),
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
    const TABLE: FileFormat = FileFormat::Table;
    const CHART: FileFormat = FileFormat::Chart;
    const FIXED: FileFormat = FileFormat::Fixed;
    match self {
      Self::All                => OutputFacts::new("",                          "--output-selection=all",       FIXED, false, "Every output of the command"),
      Self::Nwk                => OutputFacts::new(".nwk",                      "--output-tree-nwk",            FIXED, true,  "Tree in Newick format"),
      Self::Nexus              => OutputFacts::new(".nexus",                    "--output-tree-nexus",          FIXED, true,  "Tree with node annotations in Nexus format"),
      Self::Auspice            => OutputFacts::new(".auspice.json",             "--output-tree-auspice",        FIXED, true,  "Tree for Auspice and Nextstrain"),
      Self::MatPb              => OutputFacts::new(".mat.pb",                   "--output-tree-mat-pb",         FIXED, true,  "Mutation-annotated tree (UShER protobuf)"),
      Self::MatJson            => OutputFacts::new(".mat.json",                 "--output-tree-mat-json",       FIXED, true,  "Mutation-annotated tree (JSON)"),
      Self::GraphJson          => OutputFacts::new(".graph.json",               "--output-tree-graph-json",     FIXED, false, "Tree as a graph (JSON)"),
      Self::Dot                => OutputFacts::new(".dot",                      "--output-tree-dot",            FIXED, false, "Tree as a graph (Graphviz)"),
      Self::NodeData           => OutputFacts::new(".node-data.json",           "--output-node-data",           FIXED, true,  "Node data for augur: dates, intervals, mutations"),
      Self::Model              => OutputFacts::new(".model.json",               "--output-model",               FIXED, false, "Fitted substitution models"),
      Self::ClockModel         => OutputFacts::new(".clock-model.json",         "--output-clock-model",         FIXED, false, "Clock rate, intercept and regression statistics"),
      Self::Dates              => OutputFacts::new(".dates.tsv",                "--output-dates",               TABLE, true,  "Date and numeric date of every node, with the 90% interval with --confidence"),
      Self::TraitProbabilities => OutputFacts::new(".trait-probabilities.tsv",  "--output-trait-probabilities", TABLE, true,  "State probabilities of every node"),
      Self::Sequences          => OutputFacts::new(".sequences.fasta",          "--output-sequences",           FIXED, false, "Sequences of samples and ancestors"),
      Self::Translations       => OutputFacts::new(".translations.{cds}.fasta", "--output-translations",        FIXED, false, "Protein sequences of samples and ancestors"),
      Self::Traits             => OutputFacts::new(".traits.tsv",               "--output-traits",              TABLE, true,  "Inferred state of every node"),
      Self::Rtt                => OutputFacts::new(".rtt.tsv",                  "--output-rtt",                 TABLE, false, "Root-to-tip distance, date and clock prediction of every node"),
      Self::Tracelog           => OutputFacts::new(".tracelog.tsv",             "--output-tracelog",            TABLE, false, "Convergence values of every iteration"),
      Self::Coalescent         => OutputFacts::new(".coalescent.tsv",           "--output-coalescent",          TABLE, false, "Coalescent time scale and effective population size"),
      Self::CoalescentJson     => OutputFacts::new(".coalescent.json",          "--output-coalescent-json",     FIXED, false, "Coalescent model and its likelihood"),
      Self::RttChart           => OutputFacts::new(".rtt.svg",                  "--output-rtt-chart",           CHART, false, "Root-to-tip regression chart (SVG, PNG or PDF, by the file extension)"),
      Self::HomoplasyStats     => OutputFacts::new(".stats.json",               "--output-homoplasy-stats",     FIXED, false, "Recurrent mutations: counts, histograms and ranked lists (JSON)"),
      Self::HomoplasyReport    => OutputFacts::new(".report.txt",               "--output-homoplasy-report",    FIXED, false, "Recurrent mutations: text report"),
    }
  }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum TreeWriteKind {
  Nwk(NwkStyle),
  Nexus,
  Auspice,
  MatPb,
  MatJson,
  GraphJson,
  Dot,
}

impl TreeWriteKind {
  pub fn is_per_alignment(self) -> bool {
    match self {
      Self::Nwk(annotations) => annotations != NwkStyle::Plain,
      Self::Nexus | Self::Auspice | Self::MatPb | Self::MatJson => true,
      Self::GraphJson | Self::Dot => false,
    }
  }
}

struct OutputFacts {
  extension: &'static str,
  flag_name: &'static str,
  format: FileFormat,
  reads_tree: bool,
  description: &'static str,
}

impl OutputFacts {
  const fn new(
    extension: &'static str,
    flag_name: &'static str,
    format: FileFormat,
    reads_tree: bool,
    description: &'static str,
  ) -> Self {
    Self {
      extension,
      flag_name,
      format,
      reads_tree,
      description,
    }
  }
}
