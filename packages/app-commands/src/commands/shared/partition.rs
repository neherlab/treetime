use crate::commands::shared::alphabet::AlphabetNameCli;
use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::dependent_settings::check_applies;
use crate::commands::shared::documented_default::documented_default;
use crate::commands::shared::input_warnings::warn_duplicate_names;
use crate::commands::shared::model::{GtrModelNameCli, MODEL_DEFAULT, ModelSource, parse_model_params, resolve_model};
#[cfg(feature = "clap")]
use crate::yaml::yaml_read_str;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::izip;
use schemars::JsonSchema;
use serde_json::{Map, Value, json};
use treetime::alphabet::alphabet::AlphabetName;
use treetime::gtr::get_gtr::GtrModel;
use treetime::partition::resolve::{
  PartitionEntry, PartitionError, PartitionErrorKind, PartitionLayout, PartitionPlan, ResolvedPartitions,
  partition_label, resolve_partition_layouts, resolve_partitions,
};
use treetime::partition::sites::SiteRange;
use treetime::progress::{LogSink, RunWarningKind};
use treetime::progress_debug;
use treetime_schema::{schema_defaults, skip_serializing_optionals};
#[cfg(feature = "clap")]
use treetime_utils::error::report_to_string;
use treetime_utils::io::json::from_json_value;

pub fn partition_entries(specs: &[PartitionSpec]) -> Result<Vec<PartitionEntry>, Report> {
  izip!(specs, 1..).map(|(spec, number)| spec.entry(number)).collect()
}

pub fn check_partition_applies(
  partition: &[PartitionSpec],
  other: &str,
  applies: bool,
  rule: impl FnOnce() -> String,
) -> Result<(), Report> {
  check_applies("partition", !partition.is_empty(), other, applies, rule)
}

pub fn partition_argument_error(report: Report, has_entries: bool) -> Report {
  let Some(PartitionError { kind, message }) = report.downcast_ref::<PartitionError>() else {
    return report;
  };
  let settings: &[&str] = match (kind, has_entries) {
    (PartitionErrorKind::Layout, true) => &["partition", "alignment"],
    (PartitionErrorKind::Layout, false) => &["alignment"],
    (PartitionErrorKind::Model, true) => &["partition", "alphabet"],
    (PartitionErrorKind::Model, false) => &["model", "model_params", "alphabet"],
  };
  argument_error(settings, message.clone())
}

pub(crate) fn resolve_run_partitions(
  entries: &[PartitionEntry],
  model: &GtrModel,
  lengths: &[usize],
  alphabet: Option<AlphabetName>,
  log: &dyn LogSink,
) -> Result<Vec<PartitionPlan>, Report> {
  let ResolvedPartitions { plans, duplicate_names } = resolve_partitions(entries, lengths, model, alphabet)
    .map_err(|report| partition_argument_error(report, !entries.is_empty()))?;
  warn_duplicate_names(
    log,
    RunWarningKind::DuplicatePartitionNames,
    "More than one partition entry is named",
    "TreeTime tells these partitions apart by their numbers.",
    &duplicate_names,
  );
  Ok(plans)
}

pub(crate) fn resolve_alignment_layouts(lengths: &[usize]) -> Result<Vec<PartitionLayout>, Report> {
  let layouts = resolve_partition_layouts(&[], lengths).map_err(|report| partition_argument_error(report, false))?;
  Ok(layouts.layouts)
}

pub(crate) fn log_partitions(plans: &[PartitionPlan], log: &dyn LogSink) {
  for PartitionPlan { layout, model } in plans {
    progress_debug!(
      log,
      "{}: {}, model {}",
      layout.label(),
      layout_summary(layout),
      model.name().to_string().to_lowercase()
    );
  }
}

pub(crate) fn log_partition_layouts(layouts: &[PartitionLayout], log: &dyn LogSink) {
  for layout in layouts {
    progress_debug!(log, "{}: {}", layout.label(), layout_summary(layout));
  }
}

/// Partitions of the sites of the `--alignment` files, each with its own substitution model
#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct PartitionArgs {
  /// One partition: a set of sites of one `--alignment` file with its own substitution model; repeat the flag for each
  /// partition
  ///
  /// Without `--partition`, each `--alignment` file is one partition with the model of `--model`. An entry is a YAML
  /// or JSON map with the keys `name` (the name of the partition in outputs and messages), `alignment` (the number of
  /// its `--alignment` file, 1, 2, ...; required with two or more files), `sites` (a list of site ranges, numbered
  /// from 1: `a`, `a-b`, or `a-b\k` for every k-th site from a to b), `model` (the substitution model, default `infer`)
  /// and `model_params` (the `key=value` items of `--model-params` for that model). Each entry carries its own model,
  /// so `--model` and `--model-params` cannot be used together with `--partition`. Every alignment needs at least one
  /// entry, and every site belongs to exactly one entry. An entry without `sites` takes the sites that the other entries
  /// of its alignment leave; one entry per alignment can leave out `sites`. Write site ranges in single quotes or
  /// without quotes, because YAML reads a backslash in double quotes as an escape. Example: `--partition "{name: cp1,
  /// sites: ['101-1600\3']}" --partition "{name: cp2, sites: ['102-1600\3']}" --partition "{name: cp3, sites:
  /// ['103-1600\3'], model: hky85, model_params: ['kappa=4']}" --partition "{name: noncoding}"`.
  #[schemars(example = json!([
    {"name": "cp12", "sites": ["1-1500\\3", "2-1500\\3"]},
    {"name": "cp3", "sites": ["3-1500\\3"], "model": "hky85", "model_params": ["kappa=4"]},
  ]))]
  #[cfg_attr(
    feature = "clap",
    clap(
      long,
      value_name = "ENTRY",
      value_parser = parse_partition_spec,
      conflicts_with_all = ["model", "model_params"],
      hide_short_help = true,
      help_heading = "Input data",
    )
  )]
  pub partition: Vec<PartitionSpec>,
}

/// One partition: a set of sites of one `--alignment` file and its substitution model
#[derive(Debug, Clone, Default, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(deny_unknown_fields)]
#[deser(default, deny_unknown_fields)]
pub struct PartitionSpec {
  /// Name of the partition in outputs and messages
  pub name: Option<String>,

  /// Number of the alignment file in the `--alignment` list (1, 2, ...); may be left out when the list has one file
  #[schemars(range(min = 1))]
  pub alignment: Option<usize>,

  /// Sites of the partition, as site ranges numbered from 1: `a`, `a-b`, or `a-b\k` for every k-th site from a to b
  ///
  /// Left out, the partition takes every site of its alignment that no other entry takes.
  pub sites: Option<Vec<SiteRange>>,

  /// Substitution model of the partition
  #[schemars(transform = documented_default(&MODEL_DEFAULT))]
  pub model: Option<GtrModelNameCli>,

  /// Parameters of the model of the partition, as the `key=value` items of `--model-params`, for example `kappa=4`
  pub model_params: Option<Vec<String>>,
}

impl PartitionSpec {
  fn entry(&self, number: usize) -> Result<PartitionEntry, Report> {
    let label = partition_label(number, self.name.as_deref());
    let model = parse_model_params(
      self.model.unwrap_or(MODEL_DEFAULT),
      self.model_params.as_deref().unwrap_or_default(),
      ModelSource::Entry { label: &label },
    )?;
    Ok(PartitionEntry {
      name: self.name.clone(),
      alignment: self.alignment,
      sites: self.sites.clone(),
      model,
    })
  }
}

pub struct ConfigPartitions {
  entries: Vec<PartitionEntry>,
  model: GtrModel,
  alphabet: Option<AlphabetName>,
}

impl ConfigPartitions {
  pub fn from_config(config: &Map<String, Value>) -> Result<Self, Report> {
    let settings = PartitionSettings::KEYS
      .into_iter()
      .filter_map(|key| config.get(key).map(|value| (key.to_owned(), value.clone())))
      .collect();
    let PartitionSettings {
      partition,
      model,
      model_params,
      alphabet,
    } = from_json_value(&Value::Object(settings))?;
    Ok(Self {
      entries: partition_entries(&partition)?,
      model: resolve_model(model, &model_params)?,
      alphabet: alphabet.map(Into::into),
    })
  }

  pub fn resolve(&self, lengths: &[usize]) -> Result<Vec<PartitionPlan>, Report> {
    let ResolvedPartitions { plans, .. } = resolve_partitions(&self.entries, lengths, &self.model, self.alphabet)
      .map_err(|report| partition_argument_error(report, !self.entries.is_empty()))?;
    Ok(plans)
  }
}

#[derive(Default, Deserialize)]
#[deser(default, deny_unknown_fields)]
struct PartitionSettings {
  partition: Vec<PartitionSpec>,
  model: Option<GtrModelNameCli>,
  model_params: Vec<String>,
  alphabet: Option<AlphabetNameCli>,
}

impl PartitionSettings {
  const KEYS: [&str; 4] = ["partition", "model", "model_params", "alphabet"];
}

#[cfg(feature = "clap")]
fn parse_partition_spec(text: &str) -> Result<PartitionSpec, String> {
  yaml_read_str(text).map_err(|report| report_to_string(&report))
}

fn layout_summary(layout: &PartitionLayout) -> String {
  format!("alignment {}, {} sites", layout.alignment, layout.columns.len())
}
