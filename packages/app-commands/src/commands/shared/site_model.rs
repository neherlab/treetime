use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::dependent_settings::check_applies;
use crate::commands::shared::partition_by::PartitionScheme;
use crate::commands::shared::sequence_inputs::PartitionSource;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use treetime::gtr::gamma::NO_RATE_VARIATION;
use treetime::gtr::get_gtr::GtrModel;
use treetime::gtr::shared_params::{SharedParam, SharedParams};
use treetime::partition::split::ModelSettings;
use treetime_schema::{schema_defaults, skip_serializing_optionals};
use treetime_utils::fmt::string::counted;
use treetime_utils::o;

pub const MIN_GAMMA_CATEGORIES: usize = 2;

/// Rate variation across the sites of each partition, and the model parameters that partitions share
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
  /// alpha by maximum likelihood between 0.02 and 100, unless the partitions share `gamma-shape` (`--shared-params`).
  /// Without the setting, all sites of a partition share one rate. Runtime and memory grow about K-fold; 4 categories
  /// are the common choice. The setting applies to marginal reconstructions only: not with `--method-anc parsimony`
  /// and not with `--branch-length-mode input`.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, value_name = "K", help_heading = "Substitution model")
  )]
  #[schemars(range(min = 2))]
  pub gamma_categories: Option<usize>,

  /// Model parameters that the partitions with an inferred model share, as a comma-separated list
  ///
  /// `exchange-rates`: one set of exchange rates W for all these partitions. Each partition scales W so that its rate
  /// matrix has mean rate 1 under its own frequencies, so with different frequencies the stored W of the partitions
  /// differ by a factor. `frequencies`: one set of equilibrium frequencies. `gamma-shape`: one gamma shape alpha; needs
  /// `--gamma-categories`. `none`: every partition fits its own model; stands alone. The fit of the shared parameters
  /// uses the substitutions of all sharing partitions together. A partition with a named model keeps the parameters
  /// of its model. Every partition keeps its own relative rate (`--partition-rates`), and all partitions share the
  /// tree and the branch lengths.
  ///
  /// Default: with `--partition-by genes` or `rate-categories`, the partitions share `exchange-rates` and
  /// `frequencies`, and also `gamma-shape` with `--gamma-categories`, because genes of one genome share the mutation
  /// process and rate buckets differ only in speed. With `--partition-by codon-positions`, `--partition` entries, or
  /// one partition per `--alignment` file, they share nothing, because codon positions differ in their pattern of
  /// change. The setting applies only with two or more partitions whose model is `infer`.
  #[cfg_attr(
    feature = "clap",
    clap(
      long,
      hide_short_help = true,
      num_args = 1..,
      value_delimiter = ',',
      value_enum,
      help_heading = "Substitution model"
    )
  )]
  pub shared_params: Vec<SharedParamCli>,
}

/// Model parameter that partitions share, or `none`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "SharedParam")]
pub enum SharedParamCli {
  /// Share the exchange rates W.
  ExchangeRates,
  /// Share the equilibrium frequencies pi.
  Frequencies,
  /// Share the gamma shape alpha.
  GammaShape,
  /// Share no model parameter.
  None,
}

impl SharedParamCli {
  const fn param(self) -> Option<SharedParam> {
    match self {
      Self::ExchangeRates => Some(SharedParam::ExchangeRates),
      Self::Frequencies => Some(SharedParam::Frequencies),
      Self::GammaShape => Some(SharedParam::GammaShape),
      Self::None => None,
    }
  }
}

pub fn default_shared_params(scheme: Option<PartitionScheme>, gamma: bool) -> SharedParams {
  match scheme {
    Some(PartitionScheme::Genes | PartitionScheme::RateCategories) => SharedParams {
      exchange_rates: true,
      frequencies: true,
      gamma_shape: gamma,
    },
    Some(PartitionScheme::CodonPositions) | None => SharedParams::default(),
  }
}

pub(crate) fn model_settings(
  SiteModelArgs {
    gamma_categories,
    shared_params,
  }: &SiteModelArgs,
  partitions: &PartitionSource,
  model: &GtrModel,
  alignment_count: usize,
) -> Result<ModelSettings, Report> {
  let gamma_categories = gamma_categories_or_default(*gamma_categories)?;
  let gamma = gamma_categories > NO_RATE_VARIATION;
  let inferred = InferredPartitions::new(partitions, model, alignment_count);
  let shared = if shared_params.is_empty() {
    if inferred.several() {
      default_shared_params(partitions.scheme(), gamma)
    } else {
      SharedParams::default()
    }
  } else {
    explicit_shared_params(shared_params, gamma, inferred)?
  };
  Ok(ModelSettings {
    gamma_categories,
    shared,
  })
}

pub fn check_site_model_applies(
  SiteModelArgs {
    gamma_categories,
    shared_params,
  }: &SiteModelArgs,
  other: &str,
  applies: bool,
  rule: impl Fn() -> String,
) -> Result<(), Report> {
  check_applies("gamma_categories", gamma_categories.is_some(), other, applies, &rule)?;
  check_applies("shared_params", !shared_params.is_empty(), other, applies, &rule)
}

pub(crate) fn site_model_marginal_only_settings(
  SiteModelArgs {
    gamma_categories,
    shared_params,
  }: &SiteModelArgs,
) -> [(&'static str, bool); 2] {
  [
    ("gamma_categories", gamma_categories.is_some()),
    ("shared_params", !shared_params.is_empty()),
  ]
}

fn gamma_categories_or_default(gamma_categories: Option<usize>) -> Result<usize, Report> {
  match gamma_categories {
    Some(categories) if categories < MIN_GAMMA_CATEGORIES => Err(argument_error(
      &["gamma_categories"],
      format!(
        "--gamma-categories needs {MIN_GAMMA_CATEGORIES} or more categories, but got {categories}: leave the setting \
         out for a model without rate variation"
      ),
    )),
    Some(categories) => Ok(categories),
    None => Ok(NO_RATE_VARIATION),
  }
}

fn explicit_shared_params(
  shared_params: &[SharedParamCli],
  gamma: bool,
  inferred: InferredPartitions,
) -> Result<SharedParams, Report> {
  if shared_params.len() > 1 && shared_params.contains(&SharedParamCli::None) {
    return Err(argument_error(
      &["shared_params"],
      format!(
        "--shared-params takes 'none' alone, but got {}",
        shared_params
          .iter()
          .map(|param| format!("'{}'", param_name(*param)))
          .join(", ")
      ),
    ));
  }
  let shared = SharedParams::from_params(shared_params.iter().filter_map(|param| param.param()));
  check_applies("shared_params", shared.gamma_shape, "gamma_categories", gamma, || {
    o!("gamma-shape needs --gamma-categories: without rate categories, the partitions have no gamma shape")
  })?;
  inferred.check_several()?;
  Ok(shared)
}

const fn param_name(param: SharedParamCli) -> &'static str {
  match param {
    SharedParamCli::ExchangeRates => "exchange-rates",
    SharedParamCli::Frequencies => "frequencies",
    SharedParamCli::GammaShape => "gamma-shape",
    SharedParamCli::None => "none",
  }
}

#[derive(Clone, Copy, Debug)]
enum InferredPartitions {
  Entries(usize),
  Files(usize),
  Scheme,
  NamedModel,
}

impl InferredPartitions {
  fn new(partitions: &PartitionSource, model: &GtrModel, alignment_count: usize) -> Self {
    match partitions {
      PartitionSource::Entries(entries) if !entries.is_empty() => {
        Self::Entries(entries.iter().filter(|entry| entry.model.is_infer()).count())
      },
      _ if !model.is_infer() => Self::NamedModel,
      PartitionSource::Scheme(_) => Self::Scheme,
      PartitionSource::Entries(_) => Self::Files(alignment_count),
    }
  }

  const fn several(self) -> bool {
    match self {
      Self::Entries(count) | Self::Files(count) => count >= 2,
      Self::Scheme => true,
      Self::NamedModel => false,
    }
  }

  fn check_several(self) -> Result<(), Report> {
    let rule = "applies only with two or more partitions whose model is infer";
    match self {
      Self::Entries(count) => check_applies("shared_params", true, "partition", self.several(), || {
        format!(
          "{rule}, but the --partition entries have {}",
          counted(count, "such partition")
        )
      }),
      Self::Files(count) => check_applies("shared_params", true, "alignment", self.several(), || {
        format!(
          "{rule}, but the run has {}: without --partition or --partition-by, each --alignment file is one partition",
          counted(count, "partition")
        )
      }),
      Self::NamedModel => check_applies("shared_params", true, "model", false, || {
        format!("{rule}, but --model names a model, whose parameters every partition keeps")
      }),
      Self::Scheme => Ok(()),
    }
  }
}
