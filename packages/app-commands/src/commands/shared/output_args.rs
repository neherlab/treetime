use crate::commands::ancestral::args::TreetimeAncestralArgsRaw;
use crate::commands::clock::args::TreetimeClockArgsRaw;
use crate::commands::homoplasy::args::TreetimeHomoplasyArgsRaw;
use crate::commands::mugration::args::TreetimeMugrationArgsRaw;
use crate::commands::optimize::args::TreetimeOptimizeArgsRaw;
use crate::commands::prune::args::TreetimePruneArgsRaw;
use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
use app_output::output_plan::{self, CommandKind, OutputPlanRequest, OutputSelection, ResolvedOutputs};
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

macro_rules! command_outputs {
  ($name:ident for $raw:ty, $kind:ident $(alignment: $alignment:ident)? { $($extra:ident: $field:ident),* $(,)? }) => {
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

    impl ResolveOutputs for $raw {
      fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report> {
        let selection = self.output_selection.iter().copied().map(OutputSelection::from).collect_vec();
        let alignment_count = 0 $(+ self.$alignment.alignment.len())?;
        self.output.resolve(
          CommandKind::$kind,
          partition_count(alignment_count),
          &selection,
          &[$((OutputSelection::$extra, self.$field.as_deref())),*],
        )
      }
    }
  };
}

pub trait ResolveOutputs {
  fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report>;
}

fn partition_count(alignment_count: usize) -> NonZeroUsize {
  NonZeroUsize::new(alignment_count).unwrap_or(NonZeroUsize::MIN)
}

command_outputs!(AncestralOutputSelection for TreetimeAncestralArgsRaw, Ancestral alignment: alignment {
  AugurNodeData: output_augur_node_data,
  Gtr: output_gtr,
  ReconstructedNucFasta: output_reconstructed_nuc_fasta,
  ReconstructedAaFasta: output_reconstructed_aa_fasta,
});
command_outputs!(TimetreeOutputSelection for TreetimeTimetreeArgsRaw, Timetree alignment: alignment {
  AugurNodeData: output_augur_node_data,
  Gtr: output_gtr,
  ReconstructedNucFasta: output_reconstructed_nuc_fasta,
  ClockModel: output_clock_model,
  ClockCsv: output_clock_csv,
  ConfidenceTsv: output_confidence_tsv,
  Tracelog: output_tracelog,
  CoalescentTsv: output_coalescent_tsv,
  CoalescentCsv: output_coalescent_csv,
  CoalescentJson: output_coalescent_json,
});
command_outputs!(ClockOutputSelection for TreetimeClockArgsRaw, Clock {
  ClockModel: output_clock_model,
  ClockCsv: output_clock_csv,
  ClockChartSvg: output_clock_chart_svg,
  ClockChartPng: output_clock_chart_png,
});
command_outputs!(MugrationOutputSelection for TreetimeMugrationArgsRaw, Mugration {
  AugurNodeData: output_augur_node_data,
  Gtr: output_gtr,
  ConfidenceCsv: output_confidence_csv,
  TraitsCsv: output_traits_csv,
});
command_outputs!(OptimizeOutputSelection for TreetimeOptimizeArgsRaw, Optimize alignment: alignment {
  AugurNodeData: output_augur_node_data,
  Gtr: output_gtr,
});
command_outputs!(PruneOutputSelection for TreetimePruneArgsRaw, Prune alignment: alignment { Gtr: output_gtr });
command_outputs!(HomoplasyOutputSelection for TreetimeHomoplasyArgsRaw, Homoplasy {
  HomoplasyStats: output_homoplasy_stats,
  HomoplasyReport: output_homoplasy_report,
});

/// Output selection shared by every tree-writing command. `--output-all`
/// chooses the directory, `--output-selection` restricts generated files, and
/// per-file flags override or add paths. Style and topology order are separate.
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
  /// Per-file flags (`--output-tree-nwk`, `--output-augur-node-data`, etc.) override or
  /// supplement the files produced by `--output-all`.
  ///
  /// An output that the run has no data for is skipped without a message, for example the
  /// substitution model of a run that fits none. A per-file flag for such an output fails.
  ///
  /// With several `--alignment` files, the outputs that carry sequence positions are written
  /// once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER
  /// MAT, and Newick and Nexus trees in an annotated style. Files are numbered 1, 2, ... in the
  /// order of the `--alignment` list. Their default names carry the number after the command
  /// name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output
  /// must contain `{partition}`, which is replaced by the number. Other tree outputs with
  /// `{partition}` in their file name are also written once per file.
  #[cfg_attr(feature = "clap", clap(long, short = 'O', value_hint = ValueHint::DirPath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
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
  #[schemars(extend("x-path" = "output"))]
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
  #[schemars(extend("x-path" = "output"))]
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
  #[schemars(extend("x-path" = "output"))]
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
  #[schemars(extend("x-path" = "output"))]
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
    selection: &[OutputSelection],
    non_tree_fields: &[(OutputSelection, Option<&Path>)],
  ) -> Result<ResolvedOutputs, Report> {
    let request = OutputPlanRequest {
      command,
      partition_count,
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

/// CLI-facing NWK/Nexus annotation style for `--output-nwk-style`.
#[derive(Copy, Clone, Debug, Eq, PartialEq, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
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

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum DivergenceUnits {
  #[default]
  MutationsPerSite,
  Mutations,
}
