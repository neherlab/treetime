use app_output::output_plan::{OutputSelection, ResolvedOutputs, output_unavailable};
use eyre::Report;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput, ModelEntry, ModelOutput, SequenceModel, TraitModel};
use treetime::gtr::gtr::GTR;
use treetime::partition::result::PartitionResult;
use treetime::progress::LogSink;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub(crate) fn write_gtr_output(
  resolved: &ResolvedOutputs,
  models: Vec<ModelEntry>,
  missing_reason: &str,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let Some(file) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) else {
    return Ok(());
  };
  if models.is_empty() {
    return output_unavailable(OutputSelection::Gtr, file, missing_reason, log);
  }
  json_write_file(&file.path, &ModelOutput { models }, JsonPretty(true))
}

pub(crate) fn partition_models(results: &[PartitionResult], clock_rate: Option<f64>) -> Vec<ModelEntry> {
  results
    .iter()
    .filter_map(|result| {
      let gtr = result.gtr.as_ref()?;
      let model = GtrOutput::builder().gtr(gtr).model_name(result.model_name).build();
      Some(ModelEntry::Sequence(SequenceModel {
        partition: result.number,
        alignment: result.path.clone(),
        sites: result.length,
        indel_rate: result.indel_rate,
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
