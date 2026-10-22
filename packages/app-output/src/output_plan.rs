use derive_more::{Display, Error};
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use maplit::btreeset;
use schemars::JsonSchema;
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};
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
        if let Some(kind) = variant.to_tree_write_kind(request.newick_annotations) {
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
    let label = tree_write_kind_selection(kind).flag_name();
    if command.is_per_partition_tree(kind)
      || (command.reads_alignments() && PathPlaceholder::Partition.is_in_file_name(&path))
    {
      let paths = expand_path_template(label, &path, PathPlaceholder::Partition, &partitions)?;
      resolved.tree_outputs_by_partition.insert(kind, paths);
    } else {
      ensure_no_partition_placeholder(label, &path)?;
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

fn ensure_no_partition_placeholder(label: &str, path: &Path) -> Result<(), Report> {
  if PathPlaceholder::Partition.is_in_file_name(path) {
    return make_error!(
      "{label} is written once per run, so its file name cannot contain {}",
      PathPlaceholder::Partition.quoted_tokens()
    );
  }
  Ok(())
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
      Self::Partition => &["{alignment}"],
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
      Self::Partition => "alignment file",
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

pub struct OutputPlanRequest {
  pub command: CommandKind,

  pub partition_count: NonZeroUsize,

  pub output_all: Option<PathBuf>,

  pub newick_annotations: NwkStyle,

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
        label: Cow::Owned(format!("{} for alignment {partition}", selection.flag_name())),
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
    TreeWriteKind::Nexus => OutputSelection::Nexus,
    TreeWriteKind::Auspice => OutputSelection::Auspice,
    TreeWriteKind::MatPb => OutputSelection::MatPb,
    TreeWriteKind::MatJson => OutputSelection::MatJson,
    TreeWriteKind::GraphJson => OutputSelection::GraphJson,
    TreeWriteKind::Dot => OutputSelection::Dot,
  }
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
      non_tree.remove(&OutputSelection::Rtt);
      btreeset![OutputSelection::Nwk, OutputSelection::Nexus, OutputSelection::Auspice]
    } else {
      btreeset![OutputSelection::Nwk, OutputSelection::Nexus]
    };
    for non_default in [
      OutputSelection::Translations,
      OutputSelection::ConfidenceTsv,
      OutputSelection::TraitProbabilities,
      OutputSelection::Tracelog,
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
        OutputSelection::NodeData,
        OutputSelection::Model,
        OutputSelection::Sequences,
        OutputSelection::Translations
      ],
      Self::Timetree => btreeset![
        OutputSelection::NodeData,
        OutputSelection::Model,
        OutputSelection::Sequences,
        OutputSelection::ClockModel,
        OutputSelection::Rtt,
        OutputSelection::ConfidenceTsv,
        OutputSelection::Tracelog,
        OutputSelection::Coalescent,
        OutputSelection::CoalescentJson
      ],
      Self::Optimize => btreeset![OutputSelection::NodeData, OutputSelection::Model],
      Self::Mugration => btreeset![
        OutputSelection::NodeData,
        OutputSelection::Model,
        OutputSelection::TraitProbabilities,
        OutputSelection::Traits
      ],
      Self::Clock => btreeset![
        OutputSelection::ClockModel,
        OutputSelection::Rtt,
        OutputSelection::RttChart
      ],
      Self::Prune => btreeset![OutputSelection::Model],
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
      OutputSelection::Sequences => true,
      OutputSelection::NodeData => self == Self::Ancestral,
      _ => false,
    }
  }

  fn is_per_partition_tree(self, kind: TreeWriteKind) -> bool {
    self.reads_alignments() && kind.is_per_partition()
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
  ConfidenceTsv,
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
    const TSV: Option<TableFormat> = Some(TableFormat::Tsv);
    match self {
      Self::All                => OutputFacts::new("",                          "--output-selection=all",       None, false, "Every output of the command"),
      Self::Nwk                => OutputFacts::new(".nwk",                      "--output-tree-nwk",            None, true,  "Tree in Newick format"),
      Self::Nexus              => OutputFacts::new(".nexus",                    "--output-tree-nexus",          None, true,  "Tree with node annotations in Nexus format"),
      Self::Auspice            => OutputFacts::new(".auspice.json",             "--output-tree-auspice",        None, true,  "Tree for Auspice and Nextstrain"),
      Self::MatPb              => OutputFacts::new(".mat.pb",                   "--output-tree-mat-pb",         None, true,  "Mutation-annotated tree (UShER protobuf)"),
      Self::MatJson            => OutputFacts::new(".mat.json",                 "--output-tree-mat-json",       None, true,  "Mutation-annotated tree (JSON)"),
      Self::GraphJson          => OutputFacts::new(".graph.json",               "--output-tree-graph-json",     None, false, "Tree as a graph (JSON)"),
      Self::Dot                => OutputFacts::new(".dot",                      "--output-tree-dot",            None, false, "Tree as a graph (Graphviz)"),
      Self::NodeData           => OutputFacts::new(".node-data.json",           "--output-node-data",           None, true,  "Node data for augur: dates, intervals, mutations"),
      Self::Model              => OutputFacts::new(".model.json",               "--output-model",               None, false, "Fitted substitution models"),
      Self::ClockModel         => OutputFacts::new(".clock-model.json",         "--output-clock-model",         None, false, "Clock rate, intercept and regression statistics"),
      Self::ConfidenceTsv      => OutputFacts::new(".confidence.tsv",           "--output-confidence-tsv",      TSV,  false, "Date intervals of every node"),
      Self::TraitProbabilities => OutputFacts::new(".trait-probabilities.tsv",  "--output-trait-probabilities", TSV,  true,  "State probabilities of every node"),
      Self::Sequences          => OutputFacts::new(".sequences.fasta",          "--output-sequences",           None, false, "Sequences of samples and ancestors"),
      Self::Translations       => OutputFacts::new(".translations.{cds}.fasta", "--output-translations",        None, false, "Protein sequences of samples and ancestors"),
      Self::Traits             => OutputFacts::new(".traits.tsv",               "--output-traits",              TSV,  true,  "Inferred state of every node"),
      Self::Rtt                => OutputFacts::new(".rtt.tsv",                  "--output-rtt",                 TSV,  false, "Root-to-tip distance, date and clock prediction of every node"),
      Self::Tracelog           => OutputFacts::new(".tracelog.tsv",             "--output-tracelog",            TSV,  false, "Convergence values of every iteration"),
      Self::Coalescent         => OutputFacts::new(".coalescent.tsv",           "--output-coalescent",          TSV,  false, "Coalescent time scale and effective population size"),
      Self::CoalescentJson     => OutputFacts::new(".coalescent.json",          "--output-coalescent-json",     None, false, "Coalescent model and its likelihood"),
      Self::RttChart           => OutputFacts::new(".rtt.svg",                  "--output-rtt-chart",           None, false, "Root-to-tip regression chart (SVG, PNG or PDF, by the file extension)"),
      Self::HomoplasyStats     => OutputFacts::new(".stats.json",               "--output-homoplasy-stats",     None, false, "Recurrent mutations: counts, histograms and ranked lists (JSON)"),
      Self::HomoplasyReport    => OutputFacts::new(".report.txt",               "--output-homoplasy-report",    None, false, "Recurrent mutations: text report"),
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
  pub fn is_per_partition(self) -> bool {
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
