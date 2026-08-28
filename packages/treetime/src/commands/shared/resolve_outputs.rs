use crate::commands::ancestral::args::TreetimeAncestralArgs;
use crate::commands::clock::args::TreetimeClockArgs;
use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::prune::args::TreetimePruneArgs;
use crate::commands::shared::output::{CommandKind, OutputSelection, ResolvedOutputs};
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use eyre::Report;
use std::path::Path;

/// Resolve a command's configured outputs into concrete file paths.
///
/// Each command maps its per-file output flags to the shared three-tier resolver in exactly one
/// place, so the command runner and the pipeline runner (which needs a step's produced paths to
/// resolve `{{ steps.x.outputs.* }}` chaining) never drift apart.
pub trait ResolveOutputs {
  /// Output taxonomy this command draws from.
  fn command_kind(&self) -> CommandKind;

  /// Concrete output paths for this command's current configuration.
  fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report>;
}

/// Convert a command's per-command `--output-selection` list into internal selections.
fn selection<S: Copy + Into<OutputSelection>>(selection: &[S]) -> Vec<OutputSelection> {
  selection.iter().copied().map(Into::into).collect()
}

impl ResolveOutputs for TreetimeAncestralArgs {
  fn command_kind(&self) -> CommandKind {
    CommandKind::Ancestral
  }

  fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report> {
    self.output.resolve(
      CommandKind::Ancestral,
      &selection(&self.output_selection),
      &[
        (OutputSelection::AugurNodeData, self.output_augur_node_data.as_deref()),
        (OutputSelection::Gtr, self.output_gtr.as_deref()),
        (
          OutputSelection::ReconstructedNucFasta,
          self.output_reconstructed_nuc_fasta.as_deref(),
        ),
        (
          OutputSelection::ReconstructedAaFasta,
          self.output_reconstructed_aa_fasta.as_deref().map(Path::new),
        ),
      ],
    )
  }
}

impl ResolveOutputs for TreetimeTimetreeArgs {
  fn command_kind(&self) -> CommandKind {
    CommandKind::Timetree
  }

  fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report> {
    self.output.resolve(
      CommandKind::Timetree,
      &selection(&self.output_selection),
      &[
        (OutputSelection::AugurNodeData, self.output_augur_node_data.as_deref()),
        (OutputSelection::Gtr, self.output_gtr.as_deref()),
        (OutputSelection::ClockModel, self.output_clock_model.as_deref()),
        (OutputSelection::ConfidenceTsv, self.output_confidence_tsv.as_deref()),
        (OutputSelection::Tracelog, self.output_tracelog.as_deref()),
        (OutputSelection::CoalescentTsv, self.output_coalescent_tsv.as_deref()),
        (OutputSelection::CoalescentCsv, self.output_coalescent_csv.as_deref()),
        (OutputSelection::CoalescentJson, self.output_coalescent_json.as_deref()),
      ],
    )
  }
}

impl ResolveOutputs for TreetimeClockArgs {
  fn command_kind(&self) -> CommandKind {
    CommandKind::Clock
  }

  fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report> {
    self.output.resolve(
      CommandKind::Clock,
      &selection(&self.output_selection),
      &[
        (OutputSelection::ClockModel, self.output_clock_model.as_deref()),
        (OutputSelection::ClockCsv, self.output_clock_csv.as_deref()),
      ],
    )
  }
}

impl ResolveOutputs for TreetimeMugrationArgs {
  fn command_kind(&self) -> CommandKind {
    CommandKind::Mugration
  }

  fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report> {
    self.output.resolve(
      CommandKind::Mugration,
      &selection(&self.output_selection),
      &[
        (OutputSelection::AugurNodeData, self.output_augur_node_data.as_deref()),
        (OutputSelection::Gtr, self.output_gtr.as_deref()),
        (OutputSelection::ConfidenceCsv, self.output_confidence_csv.as_deref()),
        (OutputSelection::TraitsCsv, self.output_traits_csv.as_deref()),
      ],
    )
  }
}

impl ResolveOutputs for TreetimeOptimizeArgs {
  fn command_kind(&self) -> CommandKind {
    CommandKind::Optimize
  }

  fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report> {
    self.output.resolve(
      CommandKind::Optimize,
      &selection(&self.output_selection),
      &[
        (OutputSelection::AugurNodeData, self.output_augur_node_data.as_deref()),
        (OutputSelection::Gtr, self.output_gtr.as_deref()),
      ],
    )
  }
}

impl ResolveOutputs for TreetimePruneArgs {
  fn command_kind(&self) -> CommandKind {
    CommandKind::Prune
  }

  fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report> {
    self.output.resolve(
      CommandKind::Prune,
      &selection(&self.output_selection),
      &[(OutputSelection::Gtr, self.output_gtr.as_deref())],
    )
  }
}
