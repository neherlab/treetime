use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use treetime::progress::LogSink;
use treetime::progress_info;
use treetime::timetree::params::BranchLengthMode;
use treetime_graph::edge::GraphEdgeKey;

pub const AUTO_MARGINAL_MAX_BRANCH_LENGTH: f64 = 0.1;

pub fn resolve_branch_length_mode(
  mode: BranchLengthModeCli,
  has_alignment: bool,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  log: &dyn LogSink,
) -> BranchLengthMode {
  let longest = longest_branch_length(branch_lengths);
  let resolved = mode.resolve(has_alignment, longest);
  if mode == BranchLengthModeCli::Auto {
    match (has_alignment, longest) {
      (false, _) => progress_info!(log, "Branch length mode auto: no alignment, using {resolved}"),
      (true, None) => progress_info!(
        log,
        "Branch length mode auto: no input branch has a length, using {resolved}"
      ),
      (true, Some(longest)) => progress_info!(
        log,
        "Branch length mode auto: the longest input branch is {longest:.3e} substitutions per site, using {resolved}"
      ),
    }
  }
  resolved
}

pub fn longest_branch_length(branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Option<f64> {
  branch_lengths.values().flatten().copied().reduce(f64::max)
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "BranchLengthMode")]
pub enum BranchLengthModeCli {
  /// `input` without an alignment or when the longest input branch is longer than 0.1 substitutions per site,
  /// `marginal` otherwise
  #[default]
  Auto,
  /// Keep the branch lengths of the input tree
  Input,
  /// Estimate the branch lengths from the alignment by marginal ancestral reconstruction
  Marginal,
}

impl BranchLengthModeCli {
  pub fn resolve(self, has_alignment: bool, longest_branch: Option<f64>) -> BranchLengthMode {
    match self {
      Self::Input => BranchLengthMode::Input,
      Self::Marginal => BranchLengthMode::Marginal,
      Self::Auto if !has_alignment => BranchLengthMode::Input,
      Self::Auto => match longest_branch {
        Some(longest) if longest > AUTO_MARGINAL_MAX_BRANCH_LENGTH => BranchLengthMode::Input,
        _ => BranchLengthMode::Marginal,
      },
    }
  }
}
