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
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use treetime_io::nwk::NwkStyle;
#[cfg(feature = "clap")]
use treetime_schema::{schema_defaults, skip_serializing_optionals};

macro_rules! command_outputs {
  ($name:ident, $args:ident { $($extra:ident),* $(,)? }) => {
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
  };
}

pub trait RawOutputSettings {
  fn output_settings(&self) -> OutputSettings;
}

impl RawOutputSettings for TreetimeTimetreeArgsRaw {
  fn output_settings(&self) -> OutputSettings {
    timetree_output_settings(
      self.branch_length_mode,
      !self.alignment.alignment.is_empty(),
      self.coalescent.is_some(),
    )
  }
}

impl RawOutputSettings for TreetimeMugrationArgsRaw {
  fn output_settings(&self) -> OutputSettings {
    mugration_output_settings(self.confidence)
  }
}

macro_rules! default_output_settings {
  ($($raw:ty),*) => {
    $(
      impl RawOutputSettings for $raw {
        fn output_settings(&self) -> OutputSettings {
          OutputSettings::default()
        }
      }
    )*
  };
}

default_output_settings!(
  TreetimeAncestralArgsRaw,
  TreetimeClockArgsRaw,
  TreetimeHomoplasyArgsRaw,
  TreetimeOptimizeArgsRaw,
  TreetimePruneArgsRaw
);

command_outputs!(
  AncestralOutputSelection,
  AncestralOutputSelectionArgs {
    NodeData,
    Model,
    Sequences,
    Translations,
  }
);
command_outputs!(
  TimetreeOutputSelection,
  TimetreeOutputSelectionArgs {
    NodeData,
    Model,
    Sequences,
    ClockModel,
    Rtt,
    Dates,
    Tracelog,
    Coalescent,
    CoalescentJson,
  }
);
command_outputs!(
  ClockOutputSelection,
  ClockOutputSelectionArgs {
    ClockModel,
    Rtt,
    RttChart,
  }
);
command_outputs!(
  MugrationOutputSelection,
  MugrationOutputSelectionArgs {
    NodeData,
    Model,
    TraitProbabilities,
    Traits,
  }
);
command_outputs!(OptimizeOutputSelection, OptimizeOutputSelectionArgs { NodeData, Model });
command_outputs!(PruneOutputSelection, PruneOutputSelectionArgs {});
command_outputs!(
  HomoplasyOutputSelection,
  HomoplasyOutputSelectionArgs {
    HomoplasyStats,
    HomoplasyReport,
  }
);

/// Output folder, Newick node comments and tree outputs, shared by every analysis command.
#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct OutputCoreArgs {
  /// Folder that receives the default outputs of the command
  ///
  /// Writes the default outputs of the command as `<dir>/<command>.<ext>`. `--output-selection` restricts the outputs
  /// that the folder receives. The per-file flags (`--output-*`) override or supplement the files of the folder.
  ///
  /// An output that the run has no data for is skipped without a message, for example the substitution model of a run
  /// that fits none. A per-file flag for such an output fails.
  ///
  /// Every output path follows these rules. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed
  /// file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`,
  /// `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs`
  /// threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats. The
  /// path `-` of a per-file flag writes the output uncompressed to standard output. Missing parent folders are
  /// created.
  ///
  /// Table outputs choose the delimiter by the extension of their path: `.tsv` writes tabs and `.csv` commas, also
  /// before a compression suffix (`dates.tsv.gz`); `-` writes tabs to standard output. Default names use `.tsv`.
  ///
  /// With several `--alignment` files, the outputs that carry sequence positions are written once per file:
  /// reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, Nexus, and Newick with
  /// `--newick-annotations` other than `plain`. Files are numbered 1, 2, ... in the order of the `--alignment` list.
  /// Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a
  /// per-file flag for such an output must contain `{alignment}`, which is replaced by the number. Other tree outputs
  /// with `{alignment}` in their file name are also written once per file.
  ///
  /// The outputs written once per CDS, such as the reconstructed amino-acid sequences, take the placeholder `{cds}`
  /// or `%GENE` in the file name of their per-file flag, which is replaced by the name of the CDS. The placeholder is
  /// required when more than one CDS is reconstructed; folder names are used as given.
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
    clap(long, hide_short_help = true, value_enum, default_value_t = NewickAnnotationsArg::default(), help_heading = "Output")
  )]
  pub newick_annotations: NewickAnnotationsArg,

  /// Path to output Newick tree file.
  ///
  /// `--newick-annotations` sets the node comments.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tree_nwk: Option<PathBuf>,

  /// Path to output Nexus tree file.
  ///
  /// The nodes carry BEAST comments (`[&mutations="A45G,G136T",date=2016.45]`).
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tree_nexus: Option<PathBuf>,

  /// Path to output Auspice v2 JSON tree file.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tree_auspice: Option<PathBuf>,

  /// Path to output UShER MAT protobuf tree file.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tree_mat_pb: Option<PathBuf>,

  /// Path to output UShER MAT JSON tree file.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tree_mat_json: Option<PathBuf>,

  /// Path to output internal graph JSON tree file.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tree_graph_json: Option<PathBuf>,

  /// Path to output Graphviz DOT tree file.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_tree_dot: Option<PathBuf>,
}

#[derive(Debug, Clone)]
pub(crate) struct OutputRequest {
  pub(crate) command: CommandKind,
  pub(crate) alignment_count: usize,
  pub(crate) settings: OutputSettings,
  pub(crate) output_all: Option<PathBuf>,
  pub(crate) newick_annotations: NewickAnnotationsArg,
  pub(crate) selection: Vec<OutputSelection>,
  pub(crate) tree_paths: TreeOutputPaths,
  pub(crate) file_paths: Vec<(OutputSelection, Option<PathBuf>)>,
}

impl OutputRequest {
  pub(crate) fn resolve(&self) -> Result<ResolvedOutputs, Report> {
    let request = OutputPlanRequest {
      command: self.command,
      alignment_count: NonZeroUsize::new(self.alignment_count).unwrap_or(NonZeroUsize::MIN),
      output_all: self.output_all.clone(),
      newick_annotations: self.newick_annotations.into(),
      selection: self.selection.clone(),
      tree_overrides: self.tree_paths.overrides(),
      non_tree_overrides: given_paths(&self.file_paths),
      settings: self.settings,
    };
    output_plan::plan(&request).map_err(missing_outputs_argument_error)
  }
}

#[derive(Debug, Clone)]
pub(crate) struct TreeOutputPaths {
  pub(crate) nwk: Option<PathBuf>,
  pub(crate) nexus: Option<PathBuf>,
  pub(crate) auspice: Option<PathBuf>,
  pub(crate) mat_pb: Option<PathBuf>,
  pub(crate) mat_json: Option<PathBuf>,
  pub(crate) graph_json: Option<PathBuf>,
  pub(crate) dot: Option<PathBuf>,
}

impl TreeOutputPaths {
  fn overrides(&self) -> BTreeMap<OutputSelection, PathBuf> {
    given_paths(&[
      (OutputSelection::Nwk, self.nwk.clone()),
      (OutputSelection::Nexus, self.nexus.clone()),
      (OutputSelection::Auspice, self.auspice.clone()),
      (OutputSelection::MatPb, self.mat_pb.clone()),
      (OutputSelection::MatJson, self.mat_json.clone()),
      (OutputSelection::GraphJson, self.graph_json.clone()),
      (OutputSelection::Dot, self.dot.clone()),
    ])
  }
}

pub(crate) fn selected_outputs<S: Copy + Into<OutputSelection>>(selection: &[S]) -> Vec<OutputSelection> {
  selection.iter().copied().map(Into::into).collect()
}

fn given_paths(paths: &[(OutputSelection, Option<PathBuf>)]) -> BTreeMap<OutputSelection, PathBuf> {
  paths
    .iter()
    .filter_map(|(selection, path)| path.clone().map(|path| (*selection, path)))
    .collect()
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
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_clock_model: Option<PathBuf>,

  /// Path to output the root-to-tip regression table.
  ///
  /// One row per sample: the date the regression used, the root-to-tip distance it regressed on, the date the clock
  /// model predicts from that distance, and whether the clock filter excluded the sample.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
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
