#[cfg(feature = "clap")]
use clap::ValueHint;
use eyre::{Report, WrapErr};
use maplit::btreeset;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Debug;
use std::path::{Path, PathBuf};
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_graph::topology_order::{TopologyOrderPreset, TopologyOrderSpec, TopologyOrderTargetAggregate};
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::NwkStyle;
use treetime_io::nwk::{EdgeFromNwk, NodeFromNwk, nwk_read_file};
use treetime_utils::{make_error, make_report};

/// Selectable output format for `--output-selection` and the resolution system.
///
/// Tree formats are universal (available on all commands). Non-tree formats are
/// per-command (see `CommandKind::available_outputs`).
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[derive(strum_macros::EnumIter)]
pub enum OutputSelection {
  All,

  // Tree: Newick family
  Nwk,
  NwkAnnotated,
  NwkNhx,

  // Tree: Nexus family
  Nexus,
  NexusAnnotated,
  NexusNhx,

  // Tree: other
  Auspice,
  Phyloxml,
  PhyloxmlJson,
  MatPb,
  MatJson,
  GraphJson,
  Dot,

  // Non-tree
  AugurNodeData,
  Gtr,
  ClockModel,
  Confidence,
  ReconstructedNucFasta,
  ReconstructedAaFasta,
  TraitsCsv,
  ClockCsv,
}

impl OutputSelection {
  pub fn is_tree(self) -> bool {
    matches!(
      self,
      Self::Nwk
        | Self::NwkAnnotated
        | Self::NwkNhx
        | Self::Nexus
        | Self::NexusAnnotated
        | Self::NexusNhx
        | Self::Auspice
        | Self::Phyloxml
        | Self::PhyloxmlJson
        | Self::MatPb
        | Self::MatJson
        | Self::GraphJson
        | Self::Dot
    )
  }

  pub fn is_meta(self) -> bool {
    matches!(self, Self::All)
  }

  pub fn extension(self) -> &'static str {
    match self {
      Self::All => "",
      Self::Nwk => ".nwk",
      Self::NwkAnnotated => ".annotated.nwk",
      Self::NwkNhx => ".nhx.nwk",
      Self::Nexus => ".nexus",
      Self::NexusAnnotated => ".annotated.nexus",
      Self::NexusNhx => ".nhx.nexus",
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
      Self::Confidence => ".confidence.tsv",
      Self::ReconstructedNucFasta => ".reconstructed-nuc.fasta",
      Self::ReconstructedAaFasta => ".reconstructed-aa.{cds}.fasta",
      Self::TraitsCsv => ".traits.csv",
      Self::ClockCsv => ".clock.csv",
    }
  }

  pub fn to_tree_write_kind(self) -> Option<TreeWriteKind> {
    match self {
      Self::Nwk => Some(TreeWriteKind::nwk(NwkStyle::Plain)),
      Self::NwkAnnotated => Some(TreeWriteKind::nwk(NwkStyle::Beast)),
      Self::NwkNhx => Some(TreeWriteKind::nwk(NwkStyle::Nhx)),
      Self::Nexus => Some(TreeWriteKind::nexus(NwkStyle::Plain)),
      Self::NexusAnnotated => Some(TreeWriteKind::nexus(NwkStyle::Beast)),
      Self::NexusNhx => Some(TreeWriteKind::nexus(NwkStyle::Nhx)),
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
      Self::NwkAnnotated => "--output-tree-nwk-annotated",
      Self::NwkNhx => "--output-tree-nwk-nhx",
      Self::Nexus => "--output-tree-nexus",
      Self::NexusAnnotated => "--output-tree-nexus-annotated",
      Self::NexusNhx => "--output-tree-nexus-nhx",
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
      Self::Confidence => "--output-confidence",
      Self::ReconstructedNucFasta => "--output-reconstructed-nuc-fasta",
      Self::ReconstructedAaFasta => "--output-reconstructed-aa-fasta",
      Self::TraitsCsv => "--output-traits-csv",
      Self::ClockCsv => "--output-clock-csv",
    }
  }

  #[allow(clippy::option_option)]
  fn tree_field(self, args: &OutputCoreArgs) -> Option<Option<&Path>> {
    match self {
      Self::Nwk => Some(args.output_tree_nwk.as_deref()),
      Self::NwkAnnotated => Some(args.output_tree_nwk_annotated.as_deref()),
      Self::NwkNhx => Some(args.output_tree_nwk_nhx.as_deref()),
      Self::Nexus => Some(args.output_tree_nexus.as_deref()),
      Self::NexusAnnotated => Some(args.output_tree_nexus_annotated.as_deref()),
      Self::NexusNhx => Some(args.output_tree_nexus_nhx.as_deref()),
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
  #[allow(clippy::enum_glob_use)]
  pub fn available_outputs(self) -> BTreeSet<OutputSelection> {
    let tree_formats = self.available_tree_outputs();
    let non_tree = self.non_tree_outputs();
    &tree_formats | &non_tree
  }

  #[allow(clippy::enum_glob_use)]
  pub fn default_outputs(self) -> BTreeSet<OutputSelection> {
    use OutputSelection::*;
    let tree_defaults = if self == Self::Timetree {
      btreeset![Nwk, Nexus, Auspice]
    } else {
      btreeset![Nwk, Nexus]
    };
    let mut non_tree = self.non_tree_outputs();
    non_tree.remove(&ReconstructedAaFasta);
    &tree_defaults | &non_tree
  }

  #[allow(clippy::enum_glob_use)]
  fn available_tree_outputs(self) -> BTreeSet<OutputSelection> {
    use OutputSelection::*;
    let mut tree = btreeset![
      Nwk,
      NwkAnnotated,
      NwkNhx,
      Nexus,
      NexusAnnotated,
      NexusNhx,
      GraphJson,
      Dot
    ];
    if self == Self::Timetree {
      tree.insert(Auspice);
    }
    tree
  }

  #[allow(clippy::enum_glob_use)]
  fn non_tree_outputs(self) -> BTreeSet<OutputSelection> {
    use OutputSelection::*;
    match self {
      Self::Ancestral => btreeset![AugurNodeData, Gtr, ReconstructedNucFasta, ReconstructedAaFasta],
      Self::Timetree => btreeset![AugurNodeData, Gtr, ClockModel, Confidence],
      Self::Optimize => btreeset![AugurNodeData, Gtr],
      Self::Mugration => btreeset![AugurNodeData, Gtr, Confidence, TraitsCsv],
      Self::Clock => btreeset![ClockModel, ClockCsv],
      Self::Prune => btreeset![Gtr],
    }
  }

  pub fn stem(self) -> &'static str {
    match self {
      Self::Ancestral => "annotated_tree",
      Self::Timetree => "timetree",
      Self::Optimize => "annotated_tree",
      Self::Mugration => "annotated_tree",
      Self::Clock => "rerooted",
      Self::Prune => "pruned_tree",
    }
  }
}

/// Three-tier output selection shared by every tree-writing command.
///
/// Tier 1: `--output-all` bulk directory with default file names.
/// Tier 2: `--output-selection` restricts which outputs tier 1 produces.
/// Tier 3: Per-file `--output-tree-*` flags override or supplement tiers 1-2.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct OutputCoreArgs {
  /// Write all default output files into this directory.
  ///
  /// Produces the default set of tree and non-tree outputs for the command, using
  /// `<dir>/<stem>.<ext>` paths. Combine with `--output-selection` to restrict which
  /// outputs are written.
  ///
  /// Per-file flags (`--output-tree-nwk`, `--output-augur-node-data`, etc.) override or
  /// supplement the files produced by `--output-all`.
  #[cfg_attr(feature = "clap", clap(long, short = 'O', value_hint = ValueHint::DirPath))]
  pub output_all: Option<PathBuf>,

  /// Comma-separated list of outputs to produce with `--output-all`.
  ///
  /// Restricts which outputs `--output-all` writes. Special value `all` expands to the
  /// full default set. Requires `--output-all`.
  ///
  /// Per-file flags are always honored regardless of this selection.
  #[cfg_attr(feature = "clap", clap(long, value_delimiter = ',', requires = "output_all"))]
  pub output_selection: Vec<OutputSelection>,

  // -- Newick family (6 fields) --
  /// Path to output Newick tree file (plain, no annotations).
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub output_tree_nwk: Option<PathBuf>,

  /// Path to output Newick tree file with BEAST-style annotations.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub output_tree_nwk_annotated: Option<PathBuf>,

  /// Path to output Newick tree file with NHX annotations.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub output_tree_nwk_nhx: Option<PathBuf>,

  // -- Nexus family (6 fields) --
  /// Path to output Nexus tree file (plain embedded Newick).
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub output_tree_nexus: Option<PathBuf>,

  /// Path to output Nexus tree file with BEAST-style annotations.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub output_tree_nexus_annotated: Option<PathBuf>,

  /// Path to output Nexus tree file with NHX annotations.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub output_tree_nexus_nhx: Option<PathBuf>,

  // -- Other tree formats (7 fields) --
  /// Path to output Auspice v2 JSON tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub output_tree_auspice: Option<PathBuf>,

  /// Path to output PhyloXML tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub output_tree_phyloxml: Option<PathBuf>,

  /// Path to output PhyloXML-JSON tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub output_tree_phyloxml_json: Option<PathBuf>,

  /// Path to output UShER MAT protobuf tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub output_tree_mat_pb: Option<PathBuf>,

  /// Path to output UShER MAT JSON tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub output_tree_mat_json: Option<PathBuf>,

  /// Path to output internal graph JSON tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub output_tree_graph_json: Option<PathBuf>,

  /// Path to output Graphviz DOT tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub output_tree_dot: Option<PathBuf>,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub topology_order_args: TopologyOrderArgs,
}

pub struct ResolvedOutputs {
  pub tree_outputs: BTreeMap<TreeWriteKind, PathBuf>,
  pub non_tree_outputs: BTreeMap<OutputSelection, PathBuf>,
  pub topology_order: TopologyOrderSpec,
}

impl OutputCoreArgs {
  pub fn resolve<N, E, D>(
    &self,
    command: CommandKind,
    graph: &Graph<N, E, D>,
    non_tree_fields: &[(OutputSelection, Option<&Path>)],
    input_order: Option<Vec<String>>,
  ) -> Result<ResolvedOutputs, Report>
  where
    N: GraphNode + Named,
    E: GraphEdge,
    D: Sync + Send,
  {
    use strum::IntoEnumIterator;

    let stem = command.stem();
    let available = command.available_outputs();

    if !self.output_selection.is_empty() && self.output_all.is_none() {
      return make_error!("--output-selection requires --output-all");
    }

    let mut tree_outputs: BTreeMap<TreeWriteKind, PathBuf> = BTreeMap::new();
    let mut non_tree_outputs: BTreeMap<OutputSelection, PathBuf> = BTreeMap::new();

    // Collect per-file overrides: tree fields from self, non-tree from the caller
    let mut per_file_tree: BTreeMap<OutputSelection, PathBuf> = BTreeMap::new();
    let mut per_file_non_tree: BTreeMap<OutputSelection, PathBuf> = BTreeMap::new();

    for variant in OutputSelection::iter() {
      if variant.is_meta() {
        continue;
      }
      if variant.is_tree() {
        if let Some(Some(path)) = variant.tree_field(self) {
          per_file_tree.insert(variant, path.to_path_buf());
        }
      }
    }

    for &(sel, ref_path) in non_tree_fields {
      if let Some(path) = ref_path {
        per_file_non_tree.insert(sel, path.to_path_buf());
      }
    }

    // Tier 1+2: if output_all is set, compute default paths for effective selection
    if let Some(dir) = &self.output_all {
      let effective_selection =
        if self.output_selection.is_empty() || self.output_selection.contains(&OutputSelection::All) {
          command.default_outputs()
        } else {
          self.output_selection.iter().copied().collect()
        };

      for variant in &effective_selection {
        if variant.is_meta() {
          continue;
        }

        if variant.is_tree() {
          if !per_file_tree.contains_key(variant) {
            let path = dir.join(format!("{stem}{}", variant.extension()));
            per_file_tree.insert(*variant, path);
          }
        } else if !per_file_non_tree.contains_key(variant) {
          let path = dir.join(format!("{stem}{}", variant.extension()));
          per_file_non_tree.insert(*variant, path);
        }
      }
    }

    // Validate tree outputs against command availability
    for sel in per_file_tree.keys() {
      if !available.contains(sel) {
        return make_error!(
          "Output format '{}' is not available for this command. \
           Available tree outputs: {}",
          sel.flag_name(),
          available
            .iter()
            .filter(|s| s.is_tree())
            .map(|s| s.flag_name())
            .collect::<Vec<_>>()
            .join(", ")
        );
      }
    }

    // Validate non-tree outputs against command availability
    for sel in per_file_non_tree.keys() {
      if !available.contains(sel) {
        return make_error!(
          "Output format '{}' is not available for this command. \
           Available non-tree outputs: {}",
          sel.flag_name(),
          available
            .iter()
            .filter(|s| !s.is_tree() && !s.is_meta())
            .map(|s| s.flag_name())
            .collect::<Vec<_>>()
            .join(", ")
        );
      }
    }

    // Convert tree selections to TreeWriteKind
    for (sel, path) in &per_file_tree {
      if let Some(kind) = sel.to_tree_write_kind() {
        tree_outputs.insert(kind, path.clone());
      }
    }

    // Collect non-tree outputs
    for (sel, path) in per_file_non_tree {
      non_tree_outputs.insert(sel, path);
    }

    // Error if no outputs at all
    if tree_outputs.is_empty() && non_tree_outputs.is_empty() {
      return make_error!(
        "No output flags provided. At least one is required: \
         --output-all or one of the --output-tree-* / --output-* flags"
      );
    }

    // Create directories
    if let Some(dir) = &self.output_all {
      std::fs::create_dir_all(dir)
        .wrap_err_with(|| format!("Failed to create output directory '{}'", dir.display()))?;
    }
    for path in tree_outputs.values().chain(non_tree_outputs.values()) {
      if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
          std::fs::create_dir_all(parent)
            .wrap_err_with(|| format!("Failed to create parent directory '{}'", parent.display()))?;
        }
      }
    }

    // Resolve topology order
    let topology_order = self.topology_order_args.resolve_topology_order(graph, input_order)?;

    Ok(ResolvedOutputs {
      tree_outputs,
      non_tree_outputs,
      topology_order,
    })
  }
}

#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct TopologyOrderArgs {
  /// Order tree topology before writing output files.
  #[cfg_attr(feature = "clap", clap(long, value_enum))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order_target_source"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order_target_file"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order_target_aggregate"))]
  pub ladderize: Option<LadderizeArg>,

  /// Canonical topology ordering preset.
  #[cfg_attr(feature = "clap", clap(long, value_enum))]
  pub topology_order: Option<TopologyOrderArg>,

  /// Source for target-order topology sorting.
  #[cfg_attr(feature = "clap", clap(long, value_enum))]
  pub topology_order_target_source: Option<TopologyOrderTargetSourceArg>,

  /// File used by list or reference-topology target-order sources.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub topology_order_target_file: Option<PathBuf>,

  /// Aggregate used to map a subtree to a target-order position.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = TopologyOrderTargetAggregateArg::default()))]
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
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
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
