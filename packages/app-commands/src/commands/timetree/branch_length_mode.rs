use crate::commands::shared::argument_error::argument_error;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use treetime::progress::LogSink;
use treetime::progress_info;
use treetime::timetree::params::BranchLengthMode;
use treetime_graph::edge::GraphEdgeKey;
use treetime_utils::o;

pub const AUTO_MARGINAL_MAX_BRANCH_LENGTH: f64 = 0.1;

pub fn resolve_branch_length_mode(
  mode: BranchLengthModeCli,
  has_alignment: bool,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  log: &dyn LogSink,
) -> ResolvedBranchLengthMode {
  let resolved = ResolvedBranchLengthMode::new(mode, has_alignment, longest_branch_length(branch_lengths));
  if mode == BranchLengthModeCli::Auto {
    progress_info!(
      log,
      "Branch length mode auto: {}, using {}",
      resolved.reason(),
      resolved.mode
    );
  }
  resolved
}

pub fn longest_branch_length(branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Option<f64> {
  branch_lengths.values().flatten().copied().reduce(f64::max)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResolvedBranchLengthMode {
  pub given: BranchLengthModeCli,
  pub mode: BranchLengthMode,
  pub has_alignment: bool,
  pub longest_branch: Option<f64>,
}

impl ResolvedBranchLengthMode {
  pub fn new(given: BranchLengthModeCli, has_alignment: bool, longest_branch: Option<f64>) -> Self {
    Self {
      given,
      mode: given.resolve(has_alignment, longest_branch),
      has_alignment,
      longest_branch,
    }
  }

  pub fn require_marginal(&self, settings: &[&str], subject: &str) -> Result<(), Report> {
    if self.mode == BranchLengthMode::Marginal {
      return Ok(());
    }
    let settings = settings.iter().copied().chain(["branch_length_mode"]).collect_vec();
    Err(argument_error(
      &settings,
      format!(
        "{subject} requires re-estimated branch lengths (--branch-length-mode marginal), but {}",
        self.choice()
      ),
    ))
  }

  fn choice(&self) -> String {
    match self.given {
      BranchLengthModeCli::Auto => format!("--branch-length-mode auto chose {}: {}", self.mode, self.reason()),
      BranchLengthModeCli::Input | BranchLengthModeCli::Marginal => {
        format!("--branch-length-mode is {}", self.mode)
      },
    }
  }

  fn reason(&self) -> String {
    match (self.has_alignment, self.longest_branch) {
      (false, _) => o!("no alignment"),
      (true, None) => o!("no input branch has a length"),
      (true, Some(longest)) => format!("the longest input branch is {longest:.3e} substitutions per site"),
    }
  }
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
