use crate::commands::ancestral::args::TreetimeAncestralArgsRaw;
use crate::commands::clock::args::TreetimeClockArgsRaw;
use crate::commands::homoplasy::args::TreetimeHomoplasyArgsRaw;
use crate::commands::mugration::args::{TreetimeMugrationArgsRaw, mugration_output_settings};
use crate::commands::optimize::args::TreetimeOptimizeArgsRaw;
use crate::commands::prune::args::TreetimePruneArgsRaw;
use crate::commands::shared::argument_error::argument_error;
use crate::commands::timetree::args::{TreetimeTimetreeArgsRaw, timetree_output_settings};
use app_output::command_outputs::CommandKind;
use app_output::output_conditions::{OutputSettings, UnavailableOutput};
use app_output::output_format::UnsupportedExtension;
use app_output::output_plan::{
  self, AddedOutputWithoutFolder, MissingOutputs, OutputPlanRequest, OutputSelection, ResolvedOutputs,
};
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use treetime_io::nwk::NwkStyle;
#[cfg(feature = "clap")]
use treetime_schema::{schema_defaults, skip_serializing_optionals};

macro_rules! output_settings {
  ($raw:expr) => {
    OutputSettings::default()
  };
  ($raw:expr, $settings:path) => {
    $settings($raw)
  };
}

macro_rules! command_outputs {
  ($name:ident, $args:ident for $raw:ty, $kind:ident $(alignment: $alignment:ident)? $(settings: $settings:path)? { $($extra:ident: $($field:ident).+),* $(,)? }) => {
    #[derive(Copy, Clone, Debug, Eq, PartialEq, schemars::JsonSchema, strum_macros::EnumIter, Serialize, Deserialize)]
    #[schemars(rename_all = "kebab-case")]
    #[deser(rename_all = "kebab-case")]
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

    #[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
    #[schemars(default, deny_unknown_fields)]
    #[schemars(transform = schema_defaults::<Self>)]
    #[deser(default, deny_unknown_fields)]
    #[cfg_attr(feature = "clap", derive(clap::Args))]
    pub struct $args {
      /// List of outputs to produce with `--output-all`.
      ///
      /// Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this
      /// command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.
      ///
      /// A selected output that the run has no data for is skipped without a message, for example the substitution
      /// model of a run that fits none. A per-file flag for such an output fails. An output that the settings rule
      /// out fails in both forms, for example the sequences of `timetree --branch-length-mode input`.
      #[cfg_attr(
        feature = "clap",
        clap(long, short = 's', num_args = 1.., value_delimiter = ',', requires = "output_all", help_heading = "Output")
      )]
      pub output_selection: Vec<$name>,
    }

    impl ResolveOutputs for $raw {
      fn output_settings(&self) -> OutputSettings {
        output_settings!(self $(, $settings)?)
      }

      fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report> {
        let selection = self
          .output_selection
          .output_selection
          .iter()
          .copied()
          .map(OutputSelection::from)
          .collect_vec();
        let alignment_count = 0 $(+ self.$alignment.alignment.len())?;
        self.output.resolve(
          CommandKind::$kind,
          partition_count(alignment_count),
          self.output_settings(),
          &selection,
          &[$((OutputSelection::$extra, self.$($field).+.as_deref())),*],
        )
      }
    }
  };
}

pub trait ResolveOutputs {
  fn output_settings(&self) -> OutputSettings;

  fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report>;
}

fn partition_count(alignment_count: usize) -> NonZeroUsize {
  NonZeroUsize::new(alignment_count).unwrap_or(NonZeroUsize::MIN)
}

command_outputs!(AncestralOutputSelection, AncestralOutputSelectionArgs for TreetimeAncestralArgsRaw, Ancestral alignment: alignment {
  NodeData: output_node_data,
  Model: gtr_output.output_model,
  Sequences: reconstructed_fasta_output.output_sequences,
  Translations: output_translations,
});
command_outputs!(TimetreeOutputSelection, TimetreeOutputSelectionArgs for TreetimeTimetreeArgsRaw, Timetree alignment: alignment settings: timetree_output_settings {
  NodeData: output_node_data,
  Model: gtr_output.output_model,
  Sequences: reconstructed_fasta_output.output_sequences,
  ClockModel: clock_output.output_clock_model,
  Rtt: clock_output.output_rtt,
  ConfidenceTsv: output_confidence_tsv,
  Tracelog: output_tracelog,
  Coalescent: output_coalescent,
  CoalescentJson: output_coalescent_json,
});
command_outputs!(ClockOutputSelection, ClockOutputSelectionArgs for TreetimeClockArgsRaw, Clock {
  ClockModel: clock_output.output_clock_model,
  Rtt: clock_output.output_rtt,
  RttChart: output_rtt_chart,
});
command_outputs!(MugrationOutputSelection, MugrationOutputSelectionArgs for TreetimeMugrationArgsRaw, Mugration settings: mugration_output_settings {
  NodeData: output_node_data,
  Model: gtr_output.output_model,
  TraitProbabilities: output_trait_probabilities,
  Traits: output_traits,
});
command_outputs!(OptimizeOutputSelection, OptimizeOutputSelectionArgs for TreetimeOptimizeArgsRaw, Optimize alignment: alignment {
  NodeData: output_node_data,
  Model: gtr_output.output_model,
});
command_outputs!(PruneOutputSelection, PruneOutputSelectionArgs for TreetimePruneArgsRaw, Prune alignment: alignment {});
command_outputs!(HomoplasyOutputSelection, HomoplasyOutputSelectionArgs for TreetimeHomoplasyArgsRaw, Homoplasy {
  HomoplasyStats: output_homoplasy_stats,
  HomoplasyReport: output_homoplasy_report,
});

/// Output selection shared by every tree-writing command. `--output-all`
/// chooses the directory, `--output-selection` restricts generated files, and
/// per-file flags override or add paths. Newick annotations and topology order are separate.
#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct OutputCoreArgs {
  /// Write all default output files into this directory.
  ///
  /// Produces the default set of tree and non-tree outputs for the command, using
  /// `<dir>/<command>.<ext>` paths. Combine with `--output-selection` to restrict which
  /// outputs are written.
  ///
  /// Per-file flags (`--output-tree-nwk`, `--output-node-data`, etc.) override or
  /// supplement the files produced by `--output-all`.
  ///
  /// An output that the run has no data for is skipped without a message, for example the
  /// substitution model of a run that fits none. A per-file flag for such an output fails.
  ///
  /// Table outputs choose the delimiter by the extension of their path: `.tsv` writes tabs and
  /// `.csv` commas, also before a compression suffix (`dates.tsv.gz`); `-` writes tabs to
  /// standard output. Default names use `.tsv`.
  ///
  /// With several `--alignment` files, the outputs that carry sequence positions are written
  /// once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER
  /// MAT, Nexus, and Newick with `--newick-annotations` other than `plain`. Files are numbered 1, 2, ... in the
  /// order of the `--alignment` list. Their default names carry the number after the command
  /// name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output
  /// must contain `{alignment}`, which is replaced by the number. Other tree outputs with
  /// `{alignment}` in their file name are also written once per file.
  ///
  /// An output path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The
  /// compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`,
  /// `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress
  /// with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and
  /// to read than the other formats.
  #[cfg_attr(feature = "clap", clap(long, short = 'O', value_hint = ValueHint::DirPath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_all: Option<PathBuf>,

  /// Node comments of the Newick outputs: `plain`, `beast`, or `nhx`.
  ///
  /// `plain` writes no comments. `beast` writes a comment such as
  /// `[&mutations="A45G,G136T",date=2016.45]` after each node, and `nhx` writes the same values as
  /// `[&&NHX:key=value]`. The comments hold the mutations on the branch above the node, its date,
  /// and its inferred trait value, as far as the command infers them. Nexus outputs always carry
  /// `beast` comments, which FigTree reads.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_enum, default_value_t = NewickAnnotationsArg::default(), help_heading = "Output")
  )]
  pub newick_annotations: NewickAnnotationsArg,

  /// Path to output Newick tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  /// `--newick-annotations` sets the node comments.
  ///
  /// With several `--alignment` files and `--newick-annotations` other than `plain`, one file is
  /// written per alignment, and the file name must contain `{alignment}`, which is replaced by the
  /// number of the alignment (1, 2, ...).
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tree_nwk: Option<PathBuf>,

  /// Path to output Nexus tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`. The
  /// nodes carry BEAST comments (`[&mutations="A45G,G136T",date=2016.45]`).
  ///
  /// With several `--alignment` files, one file is written per alignment, and the file name must
  /// contain `{alignment}`, which is replaced by the number of the alignment (1, 2, ...).
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tree_nexus: Option<PathBuf>,

  /// Path to output Auspice v2 JSON tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// With several `--alignment` files, one file is written per alignment, and the file name
  /// must contain `{alignment}`, which is replaced by the number of the alignment (1, 2, ...).
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tree_auspice: Option<PathBuf>,

  /// Path to output UShER MAT protobuf tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// With several `--alignment` files, one file is written per alignment, and the file name
  /// must contain `{alignment}`, which is replaced by the number of the alignment (1, 2, ...).
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tree_mat_pb: Option<PathBuf>,

  /// Path to output UShER MAT JSON tree file.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  ///
  /// With several `--alignment` files, one file is written per alignment, and the file name
  /// must contain `{alignment}`, which is replaced by the number of the alignment (1, 2, ...).
  ///
  /// Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output.
  /// Use `-` to write uncompressed to stdout.
  ///
  /// Parent directories are created if missing.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
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
  #[schemars(extend("x-path" = "output"))]
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
  #[schemars(extend("x-path" = "output"))]
  pub output_tree_dot: Option<PathBuf>,
}

impl OutputCoreArgs {
  pub fn resolve(
    &self,
    command: CommandKind,
    partition_count: NonZeroUsize,
    settings: OutputSettings,
    selection: &[OutputSelection],
    non_tree_fields: &[(OutputSelection, Option<&Path>)],
  ) -> Result<ResolvedOutputs, Report> {
    let request = OutputPlanRequest {
      command,
      partition_count,
      output_all: self.output_all.clone(),
      newick_annotations: self.newick_annotations.into(),
      selection: selection.to_vec(),
      tree_overrides: self.tree_overrides(),
      non_tree_overrides: non_tree_fields
        .iter()
        .filter_map(|&(sel, path)| path.map(|path| (sel, path.to_path_buf())))
        .collect(),
      settings,
    };
    output_plan::plan(&request).map_err(missing_outputs_argument_error)
  }

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

fn missing_outputs_argument_error(report: Report) -> Report {
  if let Some(missing) = report.downcast_ref::<MissingOutputs>() {
    return match missing {
      MissingOutputs::SelectionWithoutFolder => {
        argument_error(&["output_selection", "output_all"], missing.to_string())
      },
      MissingOutputs::NoOutputs => argument_error(&["output_all"], missing.to_string()),
    };
  }
  if let Some(unsupported) = report.downcast_ref::<UnsupportedExtension>() {
    return argument_error(&[&unsupported.selection.setting_key()], unsupported.to_string());
  }
  if let Some(unavailable) = report.downcast_ref::<UnavailableOutput>() {
    return argument_error(
      &[&unavailable.request_setting(), unavailable.cause.setting()],
      unavailable.to_string(),
    );
  }
  if let Some(added) = report.downcast_ref::<AddedOutputWithoutFolder>() {
    return argument_error(&[added.setting, &added.selection.setting_key()], added.to_string());
  }
  report
}

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ClockOutputArgs {
  /// Path to output clock model JSON.
  ///
  /// Takes precedence over paths configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_clock_model: Option<PathBuf>,

  /// Path to output the root-to-tip regression table.
  ///
  /// One row per sample: the date the regression used, the root-to-tip distance it regressed on, the date the clock
  /// model predicts from that distance, and whether the clock filter excluded the sample. Takes precedence over paths
  /// configured with `--output-all` and `--output-selection`.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_rtt: Option<PathBuf>,
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum NewickAnnotationsArg {
  #[default]
  Plain,
  Beast,
  Nhx,
}

impl From<NewickAnnotationsArg> for NwkStyle {
  fn from(value: NewickAnnotationsArg) -> Self {
    match value {
      NewickAnnotationsArg::Plain => Self::Plain,
      NewickAnnotationsArg::Beast => Self::Beast,
      NewickAnnotationsArg::Nhx => Self::Nhx,
    }
  }
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum DivergenceUnits {
  #[default]
  MutationsPerSite,
  Mutations,
}
