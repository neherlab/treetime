use app_output::output_plan::{OutputSelection, ResolvedOutputs, output_unavailable};
#[cfg(feature = "clap")]
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use std::path::PathBuf;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput, ModelEntry, ModelOutput, SequenceModel, TraitModel};
use treetime::gtr::gtr::GTR;
use treetime::partition::result::{AlignmentResult, PartitionResult};
use treetime::progress::LogSink;
use treetime_schema::{schema_defaults, skip_serializing_optionals};
use treetime_utils::io::json::{JsonPretty, json_write_file};
use treetime_utils::make_internal_error;

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct GtrOutputArgs {
  /// Path to output the fitted substitution models as JSON.
  #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_model: Option<PathBuf>,
}

pub(crate) fn write_gtr_output(
  resolved: &ResolvedOutputs,
  models: Vec<ModelEntry>,
  missing_reason: &str,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let Some(file) = resolved.non_tree_outputs.get(&OutputSelection::Model) else {
    return Ok(());
  };
  if models.is_empty() {
    return output_unavailable(OutputSelection::Model, file, missing_reason, log);
  }
  json_write_file(&file.path, &ModelOutput { models }, JsonPretty(true))
}

pub(crate) fn partition_models(
  partitions: &[PartitionResult],
  alignments: &[AlignmentResult],
  clock_rate: Option<f64>,
) -> Result<Vec<ModelEntry>, Report> {
  partitions
    .iter()
    .filter_map(|partition| Some((partition, partition.gtr.as_ref()?)))
    .map(|(partition, gtr)| {
      let Some(alignment) = alignments
        .iter()
        .find(|alignment| alignment.number == partition.alignment)
      else {
        return make_internal_error!(
          "Partition {} belongs to alignment {}, which the run did not return",
          partition.number,
          partition.alignment
        );
      };
      let model = GtrOutput::builder().gtr(gtr).model_name(partition.model_name).build();
      Ok(ModelEntry::Sequence(SequenceModel {
        partition: partition.number,
        alignment: alignment.path.clone(),
        sites: partition.sites,
        indel_rate: alignment.indel_rate,
        clock_rate: clock_rate.map(|rate| rate * model.mu()),
        model,
      }))
    })
    .collect()
}

pub(crate) fn trait_model(gtr: &GTR, attribute: &str, states: Vec<String>) -> ModelEntry {
  ModelEntry::Trait(TraitModel {
    attribute: attribute.to_owned(),
    states,
    model: GtrOutput::builder().gtr(gtr).model_name(GtrModelName::Infer).build(),
  })
}
