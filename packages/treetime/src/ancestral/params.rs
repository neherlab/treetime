use crate::gtr::get_gtr::GtrModel;
use crate::partition::marginal::sample::SampleMode;

pub struct AncestralParams {
  pub method: MethodAncestral,
  pub model: GtrModel,
  pub dense: Option<bool>,
  pub include_leaves: bool,
  pub report_ambiguous: bool,
  pub impute_missing_data: bool,
  pub seed: u64,
  pub sample_from_profile: SampleMode,
}

#[derive(Copy, Debug, Clone, Default, PartialEq, Eq, PartialOrd, Ord)]
pub enum MethodAncestral {
  #[default]
  Marginal,
  Parsimony,
}
