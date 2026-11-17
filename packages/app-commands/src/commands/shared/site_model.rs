use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::dependent_settings::check_applies;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use treetime::gtr::gamma::NO_RATE_VARIATION;
use treetime::partition::split::ModelSettings;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

pub const MIN_GAMMA_CATEGORIES: usize = 2;

/// Rate variation across the sites of each partition
#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct SiteModelArgs {
  /// Number of discrete gamma rate categories of each partition (2 or more)
  ///
  /// Sites of one partition evolve at different speeds. With this setting, every partition gets K rate categories of
  /// equal probability whose rates are the means of K equal-probability slices of a gamma distribution with mean 1 and
  /// shape alpha (Yang 1994). The likelihood of a site is the mean of its likelihoods under the K rates, and the
  /// ancestral states are weighted by the probability of each rate at the site. Each partition estimates its own
  /// alpha by maximum likelihood between 0.02 and 100. Without the setting, all sites of a partition share one rate.
  /// Runtime and memory grow about K-fold; 4 categories are the common choice. The setting applies to marginal
  /// reconstructions only: not with `--method-anc parsimony` and not with `--branch-length-mode input`.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, value_name = "K", help_heading = "Substitution model")
  )]
  #[schemars(range(min = 2))]
  pub gamma_categories: Option<usize>,
}

pub fn site_model_settings(SiteModelArgs { gamma_categories }: &SiteModelArgs) -> Result<ModelSettings, Report> {
  if let Some(categories) = *gamma_categories
    && categories < MIN_GAMMA_CATEGORIES
  {
    return Err(argument_error(
      &["gamma_categories"],
      format!(
        "--gamma-categories needs {MIN_GAMMA_CATEGORIES} or more categories, but got {categories}: leave the setting \
         out for a model without rate variation"
      ),
    ));
  }
  Ok(ModelSettings {
    gamma_categories: gamma_categories.unwrap_or(NO_RATE_VARIATION),
  })
}

pub fn check_site_model_applies(
  SiteModelArgs { gamma_categories }: &SiteModelArgs,
  other: &str,
  applies: bool,
  rule: impl FnOnce() -> String,
) -> Result<(), Report> {
  check_applies("gamma_categories", gamma_categories.is_some(), other, applies, rule)
}
