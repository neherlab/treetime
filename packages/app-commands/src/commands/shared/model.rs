use crate::commands::shared::argument_error::ArgumentError;
use crate::commands::shared::dependent_settings::check_applies;
use crate::commands::shared::documented_default::documented_default;
use crate::commands::shared::mirror_enum::impl_from_mirror;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use ndarray::Array1;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use strum_macros::{Display, EnumIter, EnumString};
use treetime::gtr::get_gtr::{
  F81Params, GtrModel, GtrModelName, HKY85Params, JC69Params, Jtt92Params, K80Params, T92Params, TN93Params,
};
use treetime::optimize::partition_rates::PartitionRates;
use treetime_schema::{schema_defaults, skip_serializing_optionals};
use treetime_utils::fmt::string::counted;
use treetime_utils::o;

pub const MODEL_DEFAULT: GtrModelNameCli = GtrModelNameCli::Infer;

pub const PARTITION_RATES_DEFAULT: PartitionRatesCli = PartitionRatesCli::Estimate;

/// Substitution model selection shared by every command that infers or applies a rate matrix.
///
/// One flag name (`--model`) serves every command; the value set includes non-GTR models (for example `jtt92`).
/// `--model-params` carries model-specific `key=value` parameters.
#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
#[cfg_attr(feature = "clap", command(mut_arg("model", |arg| documented_default(&MODEL_DEFAULT).help(arg))))]
pub struct ModelArgs {
  /// Substitution model to use
  ///
  /// `infer` infers a model from the data of each partition; each `--alignment` file is one partition. Alternatively,
  /// name the model type. If the named model takes parameters, use `--model-params` to set them. The model applies to
  /// every file. In commands with `--partition`, each entry carries its own model, so `--model` is an error next to
  /// entries.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_enum, ignore_case = true, help_heading = "Substitution model")
  )]
  #[schemars(transform = documented_default(&MODEL_DEFAULT))]
  pub model: Option<GtrModelNameCli>,

  /// Parameters for the model selected by `--model`, given as space-separated `key=value` items
  ///
  /// The parameters apply to the model of every partition; each `--alignment` file is one partition. In commands
  /// with `--partition`, each entry carries its own model and parameters, so `--model-params` is an error next to
  /// entries. Keys: `mu` (the overall substitution rate, every named model), `kappa` (the transition to transversion
  /// ratio of `k80`, `hky85` and `t92`), `kappa1` and `kappa2` (the rate ratios of `tn93`), `pi` (the equilibrium
  /// frequencies of `f81`, `hky85` and `tn93`: comma-separated numbers, one frequency per state of the alphabet in its
  /// order, for example A, C, G and T for nucleotides, normalized to sum 1; also `pis`, `Pi`, `Pis`) and `pi_GC` (the
  /// GC content of `t92`, between 0 and 1). A key left out keeps its default. Example: `--model hky85 --model-params
  /// kappa=0.2 pis=0.25,0.25,0.25,0.25`.
  #[cfg_attr(
    feature = "clap",
    clap(
      long, hide_short_help = true,
      num_args = 1..,
      help_heading = "Substitution model"
    )
  )]
  pub model_params: Vec<String>,
}

pub fn resolve_model(model: Option<GtrModelNameCli>, model_params: &[String]) -> Result<GtrModel, ArgumentError> {
  parse_model_params(model.unwrap_or(MODEL_DEFAULT), model_params, ModelSource::Flags)
}

pub fn check_model_applies(
  model: Option<GtrModelNameCli>,
  model_params: &[String],
  other: &str,
  applies: bool,
  rule: impl Fn() -> String,
) -> Result<(), Report> {
  check_applies("model", model.is_some(), other, applies, &rule)?;
  check_applies("model_params", !model_params.is_empty(), other, applies, &rule)
}

pub fn parse_model_params(
  model: GtrModelNameCli,
  items: &[String],
  source: ModelSource<'_>,
) -> Result<GtrModel, ArgumentError> {
  let params = source.params_setting();
  if let (GtrModelNameCli::Infer, Some(item)) = (model, items.first()) {
    return Err(source.model_error(format!(
      "{params} needs a named {}: {model} fits the model from the data and takes no parameters, but got '{item}'",
      source.model_setting()
    )));
  }
  let accepted = accepted_params(model);
  let mut given = BTreeMap::<ModelParam, &str>::new();
  let mut values = ModelParamValues::default();
  for item in items {
    let (key, value) = item
      .split_once('=')
      .filter(|(key, value)| !key.is_empty() && !value.is_empty() && !value.contains('='))
      .ok_or_else(|| source.params_error(format!("{params} takes key=value items, but got '{item}'")))?;
    let param = key
      .parse::<ModelParam>()
      .ok()
      .filter(|param| accepted.contains(param))
      .ok_or_else(|| {
        source.model_error(format!(
          "{params} item '{item}': the model {model} has no parameter '{key}'; its parameters are {}",
          accepted.iter().join(", ")
        ))
      })?;
    if let Some(previous) = given.insert(param, item) {
      return Err(source.params_error(format!("{params} sets {param} twice: '{previous}' and '{item}'")));
    }
    values.set(param, item, value, source)?;
  }
  Ok(values.into_model(model.into()))
}

fn accepted_params(model: GtrModelNameCli) -> &'static [ModelParam] {
  match model {
    GtrModelNameCli::Infer => &[],
    GtrModelNameCli::JC69 | GtrModelNameCli::Jtt92 => &[ModelParam::Mu],
    GtrModelNameCli::K80 => &[ModelParam::Mu, ModelParam::Kappa],
    GtrModelNameCli::F81 => &[ModelParam::Mu, ModelParam::Pi],
    GtrModelNameCli::HKY85 => &[ModelParam::Mu, ModelParam::Kappa, ModelParam::Pi],
    GtrModelNameCli::T92 => &[ModelParam::Mu, ModelParam::Kappa, ModelParam::PiGc],
    GtrModelNameCli::TN93 => &[ModelParam::Mu, ModelParam::Kappa1, ModelParam::Kappa2, ModelParam::Pi],
  }
}

#[derive(Default)]
struct ModelParamValues {
  mu: Option<f64>,
  kappa: Option<f64>,
  kappa1: Option<f64>,
  kappa2: Option<f64>,
  pi: Option<Array1<f64>>,
  pi_gc: Option<f64>,
}

impl ModelParamValues {
  fn set(&mut self, param: ModelParam, item: &str, value: &str, source: ModelSource<'_>) -> Result<(), ArgumentError> {
    match param {
      ModelParam::Mu => self.mu = Some(non_negative_number(param, item, value, source)?),
      ModelParam::Kappa => self.kappa = Some(non_negative_number(param, item, value, source)?),
      ModelParam::Kappa1 => self.kappa1 = Some(non_negative_number(param, item, value, source)?),
      ModelParam::Kappa2 => self.kappa2 = Some(non_negative_number(param, item, value, source)?),
      ModelParam::Pi => self.pi = Some(frequencies(item, value, source)?),
      ModelParam::PiGc => self.pi_gc = Some(gc_content(item, value, source)?),
    }
    Ok(())
  }

  fn into_model(self, model: GtrModelName) -> GtrModel {
    let Self {
      mu,
      kappa,
      kappa1,
      kappa2,
      pi,
      pi_gc,
    } = self;
    match GtrModel::from(model) {
      GtrModel::Infer => GtrModel::Infer,
      GtrModel::JC69(default) => GtrModel::JC69(JC69Params {
        mu: mu.unwrap_or(default.mu),
        ..default
      }),
      GtrModel::K80(default) => GtrModel::K80(K80Params {
        mu: mu.unwrap_or(default.mu),
        kappa: kappa.unwrap_or(default.kappa),
        ..default
      }),
      GtrModel::F81(default) => GtrModel::F81(F81Params {
        mu: mu.unwrap_or(default.mu),
        pi: pi.or(default.pi),
        ..default
      }),
      GtrModel::HKY85(default) => GtrModel::HKY85(HKY85Params {
        mu: mu.unwrap_or(default.mu),
        kappa: kappa.unwrap_or(default.kappa),
        pi: pi.or(default.pi),
        ..default
      }),
      GtrModel::T92(default) => GtrModel::T92(T92Params {
        mu: mu.unwrap_or(default.mu),
        kappa: kappa.unwrap_or(default.kappa),
        pi_GC: pi_gc.unwrap_or(default.pi_GC),
        ..default
      }),
      GtrModel::TN93(default) => GtrModel::TN93(TN93Params {
        mu: mu.unwrap_or(default.mu),
        kappa1: kappa1.unwrap_or(default.kappa1),
        kappa2: kappa2.unwrap_or(default.kappa2),
        pi: pi.or(default.pi),
        ..default
      }),
      GtrModel::Jtt92(default) => GtrModel::Jtt92(Jtt92Params {
        mu: mu.unwrap_or(default.mu),
        ..default
      }),
    }
  }
}

fn non_negative_number(
  param: ModelParam,
  item: &str,
  value: &str,
  source: ModelSource<'_>,
) -> Result<f64, ArgumentError> {
  parse_number(value)
    .filter(|number| number.is_finite() && *number >= 0.0)
    .ok_or_else(|| {
      source.params_error(format!(
        "{} item '{item}': {param} must be a non-negative finite number",
        source.params_setting()
      ))
    })
}

fn gc_content(item: &str, value: &str, source: ModelSource<'_>) -> Result<f64, ArgumentError> {
  parse_number(value)
    .filter(|number| *number > 0.0 && *number < 1.0)
    .ok_or_else(|| {
      source.params_error(format!(
        "{} item '{item}': {} must be a number between 0 and 1, exclusive",
        source.params_setting(),
        ModelParam::PiGc
      ))
    })
}

fn frequencies(item: &str, value: &str, source: ModelSource<'_>) -> Result<Array1<f64>, ArgumentError> {
  value
    .split(',')
    .map(|number| {
      parse_number(number)
        .filter(|number| number.is_finite() && *number > 0.0)
        .ok_or_else(|| {
          source.params_error(format!(
            "{} item '{item}': each frequency of {} must be a positive finite number",
            source.params_setting(),
            ModelParam::Pi
          ))
        })
    })
    .collect()
}

#[cfg_attr(
  dylint_lib = "custom",
  expect(
    result_defaulted,
    reason = "the parse error of a number says only that the text is not a number, which the message of every caller states"
  )
)]
fn parse_number(text: &str) -> Option<f64> {
  text.parse::<f64>().ok()
}

#[derive(Clone, Copy, Debug)]
pub enum ModelSource<'a> {
  Flags,
  Entry { label: &'a str },
}

impl ModelSource<'_> {
  fn params_setting(self) -> String {
    match self {
      Self::Flags => o!("--model-params"),
      Self::Entry { label } => format!("{label}: model_params"),
    }
  }

  const fn model_setting(self) -> &'static str {
    match self {
      Self::Flags => "--model",
      Self::Entry { .. } => "model",
    }
  }

  fn model_error(self, message: String) -> ArgumentError {
    self.error(message, &["model", "model_params"])
  }

  fn params_error(self, message: String) -> ArgumentError {
    self.error(message, &["model_params"])
  }

  fn error(self, message: String, flag_settings: &[&str]) -> ArgumentError {
    let settings = match self {
      Self::Flags => flag_settings,
      Self::Entry { .. } => &["partition"],
    };
    ArgumentError {
      message,
      settings: settings.iter().map(|setting| (*setting).to_owned()).collect(),
    }
  }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Display, EnumString)]
enum ModelParam {
  #[strum(serialize = "mu")]
  Mu,
  #[strum(serialize = "kappa")]
  Kappa,
  #[strum(serialize = "kappa1")]
  Kappa1,
  #[strum(serialize = "kappa2")]
  Kappa2,
  #[strum(to_string = "pi", serialize = "pis", serialize = "Pi", serialize = "Pis")]
  Pi,
  #[strum(serialize = "pi_GC")]
  PiGc,
}

#[derive(
  Copy,
  Debug,
  Clone,
  PartialEq,
  Eq,
  PartialOrd,
  Ord,
  SmartDefault,
  JsonSchema,
  Serialize,
  Deserialize,
  Display,
  EnumIter,
)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "GtrModelName")]
#[strum(serialize_all = "lowercase")]
pub enum GtrModelNameCli {
  /// Infer GTR parameters from data via Fitch parsimony substitution counts.
  #[default]
  Infer,
  #[schemars(rename = "jc69")]
  #[deser(rename = "jc69")]
  JC69,
  K80,
  F81,
  #[schemars(rename = "hky85")]
  #[deser(rename = "hky85")]
  HKY85,
  T92,
  #[schemars(rename = "tn93")]
  #[deser(rename = "tn93")]
  TN93,
  #[cfg_attr(feature = "clap", value(name = "jtt92"))]
  Jtt92,
}

impl_from_mirror!(GtrModelNameCli => GtrModelName { Infer, JC69, K80, F81, HKY85, T92, TN93, Jtt92 });

/// Relative rates of the partitions, shared by the commands that optimize branch lengths.
#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
#[cfg_attr(
  feature = "clap",
  command(mut_arg("partition_rates", |arg| documented_default(&PARTITION_RATES_DEFAULT).help(arg)))
)]
pub struct PartitionRatesArgs {
  /// Relative substitution rates of the partitions
  ///
  /// Each `--alignment` file is one partition, or each `--partition` entry when there are entries. With two or more
  /// partitions, all partitions share the branch lengths, and each partition has one relative rate that scales them for
  /// the sites of that partition. `estimate` fits the rate of each partition with the rest of its substitution model
  /// fixed (a partition with `--model infer` starts from its fitted rate), then divides all rates by their mean
  /// weighted by the number of sites of each partition, so that this mean is 1. The branch lengths are then expected
  /// substitutions per site, averaged over all sites of all partitions. `equal` fixes the rate of every partition at 1.
  ///
  /// The setting applies only with two or more partitions.
  #[cfg_attr(
    feature = "clap",
    clap(long, hide_short_help = true, value_enum, help_heading = "Substitution model")
  )]
  #[schemars(transform = documented_default(&PARTITION_RATES_DEFAULT))]
  pub partition_rates: Option<PartitionRatesCli>,
}

pub fn resolve_partition_rates(partition_rates: Option<PartitionRatesCli>) -> PartitionRates {
  partition_rates.unwrap_or(PARTITION_RATES_DEFAULT).into()
}

pub fn check_partition_rates(
  partition_rates: Option<PartitionRatesCli>,
  alignment_files: usize,
  partition_entries: usize,
) -> Result<(), Report> {
  if partition_entries > 0 {
    return check_applies(
      "partition_rates",
      partition_rates.is_some(),
      "partition",
      partition_entries >= 2,
      || o!("applies only with two or more partitions, but --partition has one entry, the one partition of the run"),
    );
  }
  check_applies(
    "partition_rates",
    partition_rates.is_some(),
    "alignment",
    alignment_files >= 2,
    || {
      format!(
        "applies only with two or more partitions, but the run has {}: without --partition, each --alignment file \
         is one partition",
        counted(alignment_files, "partition")
      )
    },
  )
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "PartitionRates")]
pub enum PartitionRatesCli {
  /// Estimate one relative rate per partition, normalized to a site-weighted mean of 1.
  #[default]
  Estimate,
  /// Fix the relative rate of every partition at 1.
  Equal,
}

impl_from_mirror!(PartitionRatesCli => PartitionRates { Estimate, Equal });
