use crate::gtr::get_gtr::GtrModelName;
use crate::partition::marginal::sample::SampleMode;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use strum_macros::VariantNames;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, VariantNames, JsonSchema)]
#[strum(serialize_all = "kebab-case")]
#[derive(Default)]
#[serde(rename_all = "kebab-case")]
pub enum MethodAncestral {
  #[default]
  Marginal,
  Parsimony,
  Joint,
}

pub struct AncestralParams {
  pub method: MethodAncestral,
  pub model: GtrModelName,
  pub dense: Option<bool>,
  pub include_leaves: bool,
  pub impute_missing_data: bool,
  pub gtr_iterations: usize,
  pub site_specific_gtr: bool,
  pub seed: Option<u64>,
  pub sample_from_profile: SampleMode,
}
