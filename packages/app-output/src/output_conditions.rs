use crate::output_plan::{OutputSelection, Requested};
use derive_more::{Display, Error};

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct OutputSettings {
  pub coalescent: bool,
  pub trait_probabilities: bool,
  pub sequences: SequenceOutputs,
}

impl OutputSettings {
  pub(crate) const fn unavailable(self, selection: OutputSelection) -> Option<NotReconstructed> {
    match (selection, self.sequences) {
      (OutputSelection::Sequences, SequenceOutputs::NotReconstructed(cause)) => Some(cause),
      _ => None,
    }
  }

  pub(crate) const fn condition(self, selection: OutputSelection) -> Option<OutputCondition> {
    match (selection, self.sequences) {
      (OutputSelection::Sequences, SequenceOutputs::WhenBranchLengthsReestimated) => {
        Some(OutputCondition::BranchLengthsReestimated)
      },
      _ => None,
    }
  }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SequenceOutputs {
  #[default]
  Reconstructed,
  WhenBranchLengthsReestimated,
  NotReconstructed(NotReconstructed),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NotReconstructed {
  InputBranchLengths,
  NoAlignment,
}

impl NotReconstructed {
  pub const fn setting(self) -> &'static str {
    match self {
      Self::InputBranchLengths => "branch_length_mode",
      Self::NoAlignment => "alignment",
    }
  }

  const fn reason(self) -> &'static str {
    match self {
      Self::InputBranchLengths => {
        "--branch-length-mode input keeps the input branch lengths and reconstructs no sequences"
      },
      Self::NoAlignment => "a run without --alignment reconstructs no sequences",
    }
  }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum OutputCondition {
  BranchLengthsReestimated,
}

impl OutputCondition {
  pub const fn description(self) -> &'static str {
    match self {
      Self::BranchLengthsReestimated => "written when auto resolves to marginal",
    }
  }

  pub const fn requirement(self) -> &'static str {
    match self {
      Self::BranchLengthsReestimated => "branch_length_mode: marginal",
    }
  }

  pub(crate) const fn unmet_reason(self) -> &'static str {
    match self {
      Self::BranchLengthsReestimated => {
        "--branch-length-mode auto resolved to input, which keeps the input branch lengths and reconstructs no sequences"
      },
    }
  }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Display, Error)]
#[display("{}, but {}", request_text(*selection, *requested), cause.reason())]
pub struct UnavailableOutput {
  pub selection: OutputSelection,
  pub requested: Requested,
  pub cause: NotReconstructed,
}

impl UnavailableOutput {
  pub fn request_setting(&self) -> String {
    match self.requested {
      Requested::Named => self.selection.setting_key(),
      Requested::All => "output_selection".to_owned(),
    }
  }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Display, Error)]
#[display("{} was requested, but {}", selection.flag_name(), condition.unmet_reason())]
pub struct UnmetCondition {
  pub selection: OutputSelection,
  pub condition: OutputCondition,
}

fn request_text(selection: OutputSelection, requested: Requested) -> String {
  match requested {
    Requested::Named => format!("{} was requested", selection.flag_name()),
    Requested::All => format!("--output-selection lists {}", selection.as_ref()),
  }
}
