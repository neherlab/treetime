use crate::gtr::get_gtr::GtrModelName;
use crate::partition::marginal::sample::SampleMode;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

#[derive(Copy, Debug, Clone, Default, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "kebab-case")]
pub enum MethodAncestral {
  #[default]
  Marginal,
  Parsimony,
}

pub struct AncestralParams {
  pub method: MethodAncestral,
  pub model: GtrModelName,
  pub dense: Option<bool>,
  pub include_leaves: bool,
  pub report_ambiguous: bool,
  pub impute_missing_data: bool,
  pub gtr_iterations: usize,
  pub site_specific_gtr: bool,
  pub seed: u64,
  pub sample_from_profile: SampleMode,
}
