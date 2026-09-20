use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
use crate::commands::clock::args::{TreetimeClockArgs, TreetimeClockArgsRaw};
use crate::commands::mugration::args::{TreetimeMugrationArgs, TreetimeMugrationArgsRaw};
use crate::commands::optimize::args::{TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
use crate::commands::prune::args::{TreetimePruneArgs, TreetimePruneArgsRaw};
use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
use eyre::Report;
use std::path::Path;

pub trait ResolveOutputs {
  fn command_kind(&self) -> CommandKind;

  fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report>;
}

fn selection<S: Copy + Into<OutputSelection>>(selection: &[S]) -> Vec<OutputSelection> {
  selection.iter().copied().map(Into::into).collect()
}

macro_rules! impl_resolve_outputs {
  ($kind:ident; $($ty:ty),+ $(,)?; |$s:ident| $files:expr) => {
    $(
      impl ResolveOutputs for $ty {
        fn command_kind(&self) -> CommandKind {
          CommandKind::$kind
        }

        fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report> {
          let $s = self;
          $s.output
            .resolve(CommandKind::$kind, &selection(&$s.output_selection), &$files)
        }
      }
    )+
  };
}

impl_resolve_outputs!(Ancestral; TreetimeAncestralArgs, TreetimeAncestralArgsRaw; |s| [
  (OutputSelection::AugurNodeData, s.output_augur_node_data.as_deref()),
  (OutputSelection::Gtr, s.output_gtr.as_deref()),
  (OutputSelection::ReconstructedNucFasta, s.output_reconstructed_nuc_fasta.as_deref()),
  (OutputSelection::ReconstructedAaFasta, s.output_reconstructed_aa_fasta.as_deref().map(Path::new)),
]);

impl_resolve_outputs!(Timetree; TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw; |s| [
  (OutputSelection::AugurNodeData, s.output_augur_node_data.as_deref()),
  (OutputSelection::Gtr, s.output_gtr.as_deref()),
  (OutputSelection::ReconstructedNucFasta, s.output_reconstructed_nuc_fasta.as_deref()),
  (OutputSelection::ClockModel, s.output_clock_model.as_deref()),
  (OutputSelection::ConfidenceTsv, s.output_confidence_tsv.as_deref()),
  (OutputSelection::Tracelog, s.output_tracelog.as_deref()),
  (OutputSelection::CoalescentTsv, s.output_coalescent_tsv.as_deref()),
  (OutputSelection::CoalescentCsv, s.output_coalescent_csv.as_deref()),
  (OutputSelection::CoalescentJson, s.output_coalescent_json.as_deref()),
]);

impl_resolve_outputs!(Clock; TreetimeClockArgs, TreetimeClockArgsRaw; |s| [
  (OutputSelection::ClockModel, s.output_clock_model.as_deref()),
  (OutputSelection::ClockCsv, s.output_clock_csv.as_deref()),
]);

impl_resolve_outputs!(Mugration; TreetimeMugrationArgs, TreetimeMugrationArgsRaw; |s| [
  (OutputSelection::AugurNodeData, s.output_augur_node_data.as_deref()),
  (OutputSelection::Gtr, s.output_gtr.as_deref()),
  (OutputSelection::ConfidenceCsv, s.output_confidence_csv.as_deref()),
  (OutputSelection::TraitsCsv, s.output_traits_csv.as_deref()),
]);

impl_resolve_outputs!(Optimize; TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw; |s| [
  (OutputSelection::AugurNodeData, s.output_augur_node_data.as_deref()),
  (OutputSelection::Gtr, s.output_gtr.as_deref()),
]);

impl_resolve_outputs!(Prune; TreetimePruneArgs, TreetimePruneArgsRaw; |s| [
  (OutputSelection::Gtr, s.output_gtr.as_deref()),
]);
