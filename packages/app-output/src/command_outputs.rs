use crate::output_conditions::OutputSettings;
use crate::output_plan::{OutputSelection, TreeWriteKind};
use std::collections::BTreeSet;

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
  pub fn all_selectable(self, settings: OutputSettings) -> BTreeSet<OutputSelection> {
    self
      .outputs()
      .iter()
      .map(|&(selection, _)| selection)
      .filter(|&selection| settings.unavailable(selection).is_none())
      .collect()
  }

  pub fn default_outputs(self, settings: OutputSettings) -> BTreeSet<OutputSelection> {
    self
      .outputs()
      .iter()
      .filter(|(selection, rule)| rule.includes(settings) && settings.unavailable(*selection).is_none())
      .map(|&(selection, _)| selection)
      .collect()
  }

  pub fn added_outputs(self, settings: OutputSettings) -> impl Iterator<Item = (OutputSelection, &'static str)> {
    self
      .outputs()
      .iter()
      .filter_map(move |&(selection, rule)| Some((selection, rule.added_by(settings)?)))
  }

  pub const fn outputs(self) -> &'static [(OutputSelection, DefaultRule)] {
    match self {
      Self::Ancestral => &[
        (OutputSelection::Nwk, DefaultRule::Default),
        (OutputSelection::Nexus, DefaultRule::Default),
        (OutputSelection::Auspice, DefaultRule::OnRequest),
        (OutputSelection::MatPb, DefaultRule::OnRequest),
        (OutputSelection::MatJson, DefaultRule::OnRequest),
        (OutputSelection::GraphJson, DefaultRule::OnRequest),
        (OutputSelection::Dot, DefaultRule::OnRequest),
        (OutputSelection::NodeData, DefaultRule::Default),
        (OutputSelection::Model, DefaultRule::Default),
        (OutputSelection::Sequences, DefaultRule::Default),
        (OutputSelection::Translations, DefaultRule::OnRequest),
      ],
      Self::Timetree => &[
        (OutputSelection::Nwk, DefaultRule::Default),
        (OutputSelection::Nexus, DefaultRule::Default),
        (OutputSelection::Auspice, DefaultRule::Default),
        (OutputSelection::MatPb, DefaultRule::OnRequest),
        (OutputSelection::MatJson, DefaultRule::OnRequest),
        (OutputSelection::GraphJson, DefaultRule::OnRequest),
        (OutputSelection::Dot, DefaultRule::OnRequest),
        (OutputSelection::NodeData, DefaultRule::Default),
        (OutputSelection::Model, DefaultRule::Default),
        (OutputSelection::Sequences, DefaultRule::Default),
        (OutputSelection::ClockModel, DefaultRule::Default),
        (OutputSelection::Rtt, DefaultRule::OnRequest),
        (OutputSelection::Dates, DefaultRule::Default),
        (OutputSelection::Tracelog, DefaultRule::OnRequest),
        (OutputSelection::Coalescent, DefaultRule::WithCoalescent),
        (OutputSelection::CoalescentJson, DefaultRule::OnRequest),
      ],
      Self::Optimize => &[
        (OutputSelection::Nwk, DefaultRule::Default),
        (OutputSelection::Nexus, DefaultRule::Default),
        (OutputSelection::Auspice, DefaultRule::OnRequest),
        (OutputSelection::MatPb, DefaultRule::OnRequest),
        (OutputSelection::MatJson, DefaultRule::OnRequest),
        (OutputSelection::GraphJson, DefaultRule::OnRequest),
        (OutputSelection::Dot, DefaultRule::OnRequest),
        (OutputSelection::NodeData, DefaultRule::Default),
        (OutputSelection::Model, DefaultRule::Default),
      ],
      Self::Mugration => &[
        (OutputSelection::Nwk, DefaultRule::Default),
        (OutputSelection::Nexus, DefaultRule::Default),
        (OutputSelection::Auspice, DefaultRule::OnRequest),
        (OutputSelection::MatPb, DefaultRule::OnRequest),
        (OutputSelection::MatJson, DefaultRule::OnRequest),
        (OutputSelection::GraphJson, DefaultRule::OnRequest),
        (OutputSelection::Dot, DefaultRule::OnRequest),
        (OutputSelection::NodeData, DefaultRule::Default),
        (OutputSelection::Model, DefaultRule::Default),
        (OutputSelection::Traits, DefaultRule::Default),
        (OutputSelection::TraitProbabilities, DefaultRule::WithConfidence),
      ],
      Self::Clock => &[
        (OutputSelection::Nwk, DefaultRule::Default),
        (OutputSelection::Nexus, DefaultRule::Default),
        (OutputSelection::Auspice, DefaultRule::OnRequest),
        (OutputSelection::MatPb, DefaultRule::OnRequest),
        (OutputSelection::MatJson, DefaultRule::OnRequest),
        (OutputSelection::GraphJson, DefaultRule::OnRequest),
        (OutputSelection::Dot, DefaultRule::OnRequest),
        (OutputSelection::ClockModel, DefaultRule::Default),
        (OutputSelection::Rtt, DefaultRule::Default),
        (OutputSelection::RttChart, DefaultRule::Default),
      ],
      Self::Prune => &[
        (OutputSelection::Nwk, DefaultRule::Default),
        (OutputSelection::Nexus, DefaultRule::Default),
        (OutputSelection::Auspice, DefaultRule::OnRequest),
        (OutputSelection::MatPb, DefaultRule::OnRequest),
        (OutputSelection::MatJson, DefaultRule::OnRequest),
        (OutputSelection::GraphJson, DefaultRule::OnRequest),
        (OutputSelection::Dot, DefaultRule::OnRequest),
      ],
      Self::Homoplasy => &[
        (OutputSelection::Nwk, DefaultRule::Default),
        (OutputSelection::Nexus, DefaultRule::Default),
        (OutputSelection::Auspice, DefaultRule::OnRequest),
        (OutputSelection::MatPb, DefaultRule::OnRequest),
        (OutputSelection::MatJson, DefaultRule::OnRequest),
        (OutputSelection::GraphJson, DefaultRule::OnRequest),
        (OutputSelection::Dot, DefaultRule::OnRequest),
        (OutputSelection::HomoplasyStats, DefaultRule::Default),
        (OutputSelection::HomoplasyReport, DefaultRule::Default),
      ],
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

  pub(crate) fn reads_alignments(self) -> bool {
    match self {
      Self::Ancestral | Self::Timetree | Self::Optimize | Self::Prune | Self::Homoplasy => true,
      Self::Mugration | Self::Clock => false,
    }
  }

  pub(crate) fn is_per_alignment(self, selection: OutputSelection) -> bool {
    match selection {
      OutputSelection::Sequences => true,
      OutputSelection::NodeData => self == Self::Ancestral,
      _ => false,
    }
  }

  pub(crate) fn is_per_alignment_tree(self, kind: TreeWriteKind) -> bool {
    self.reads_alignments() && kind.is_per_alignment()
  }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DefaultRule {
  Default,
  OnRequest,
  WithCoalescent,
  WithConfidence,
}

impl DefaultRule {
  const fn includes(self, settings: OutputSettings) -> bool {
    match self {
      Self::Default => true,
      Self::OnRequest => false,
      Self::WithCoalescent => settings.coalescent,
      Self::WithConfidence => settings.trait_probabilities,
    }
  }

  const fn added_by(self, settings: OutputSettings) -> Option<&'static str> {
    match self {
      Self::WithConfidence if settings.trait_probabilities => Some("confidence"),
      Self::Default | Self::OnRequest | Self::WithCoalescent | Self::WithConfidence => None,
    }
  }
}
