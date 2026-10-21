use crate::output_plan::{OutputSelection, ResolvedOutputs};
use std::collections::{BTreeMap, BTreeSet};
use treetime::partition::marginal::posterior_states::AmbiguityThreshold;
use treetime::progress::LogSink;
use treetime::progress_warn;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AmbiguityOutputs {
  carriers: BTreeMap<OutputSelection, AmbiguityThreshold>,
  selected: BTreeSet<OutputSelection>,
}

impl AmbiguityOutputs {
  pub fn new(
    threshold: AmbiguityThreshold,
    selected: &BTreeSet<OutputSelection>,
    carriers: impl Iterator<Item = OutputSelection>,
  ) -> Self {
    let carriers = carriers
      .filter(|output| selected.is_empty() || selected.contains(output))
      .map(|output| (output, threshold))
      .collect();
    Self {
      carriers,
      selected: selected.clone(),
    }
  }

  pub fn threshold_for(&self, output: OutputSelection) -> Option<AmbiguityThreshold> {
    self.carriers.get(&output).copied()
  }

  pub fn is_requested(&self) -> bool {
    !self.carriers.is_empty()
  }

  pub fn warn_unwritten(&self, resolved: &ResolvedOutputs, log: &dyn LogSink) {
    let written = resolved.paths_by_selection();
    for output in self.selected.iter().filter(|output| !written.contains_key(output)) {
      progress_warn!(
        log,
        "The ambiguity report is requested for the output '{}' ({}), but the run does not write this output",
        output.as_ref(),
        output.flag_name()
      );
    }
  }
}
