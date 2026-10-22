use crate::commands::shared::dependent_settings::check_applies;
use crate::commands::shared::mirror_enum::impl_from_mirror;
use crate::commands::shared::model::{GtrModelNameCli, check_model_applies};
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use treetime::ancestral::params::MethodAncestral;
use treetime::partition::marginal::sample::SampleMode;
use treetime_schema::schema_defaults;
use treetime_utils::o;

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct MethodAncArgs {
  /// Method used for reconstructing ancestral sequences, which places the mutations on the branches
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = MethodAncestralCli::default(), help_heading = "Ancestral reconstruction"))]
  pub method_anc: MethodAncestralCli,
}

pub fn check_model_with_method_anc(
  method_anc: MethodAncestralCli,
  model: Option<GtrModelNameCli>,
  model_params: &[String],
) -> Result<(), Report> {
  check_model_applies(model, model_params, "method_anc", method_anc.is_marginal(), || {
    o!("applies only with --method-anc marginal: parsimony reconstruction uses no substitution model")
  })
}

pub struct MarginalReconstruction {
  pub runs: bool,
  pub condition: &'static str,
}

impl MarginalReconstruction {
  pub fn check(&self, setting: &str, given: bool) -> Result<(), Report> {
    check_applies(setting, given, "method_anc", self.runs, || {
      format!(
        "applies only when a marginal reconstruction runs ({}), but --method-anc is parsimony",
        self.condition
      )
    })
  }
}

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct SampleFromProfileArgs {
  /// How to pick ancestral states from the marginal posterior profile.
  ///
  /// 'argmax': most likely state at every node (deterministic, default).
  /// 'root': sample from the posterior at the root only, argmax elsewhere (as augur's `sample_from_profile='root'`).
  /// 'all': sample from the posterior at every node.
  ///
  /// Only affects marginal reconstruction (`--method-anc=marginal`). Use `--seed` for reproducible draws.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = SampleModeCli::default(), help_heading = "Ancestral reconstruction"))]
  pub sample_from_profile: SampleModeCli,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Default, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "MethodAncestral")]
pub enum MethodAncestralCli {
  #[default]
  Marginal,
  Parsimony,
}

impl MethodAncestralCli {
  pub fn is_marginal(self) -> bool {
    self == Self::Marginal
  }
}

impl_from_mirror!(MethodAncestralCli => MethodAncestral { Marginal, Parsimony });

#[derive(Clone, Copy, Debug, PartialEq, Eq, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "SampleMode")]
pub enum SampleModeCli {
  #[default]
  Argmax,
  Root,
  All,
}

impl_from_mirror!(SampleModeCli => SampleMode { Argmax, Root, All });
