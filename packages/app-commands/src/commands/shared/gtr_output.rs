use app_output::output_plan::{OutputSelection, ResolvedOutputs, output_unavailable};
use eyre::Report;
use treetime::gtr::get_gtr::{GtrModelName, GtrOutput};
use treetime::gtr::gtr::GTR;
use treetime::progress::LogSink;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub(crate) fn write_gtr_output(
  resolved: &ResolvedOutputs,
  model: Option<(&GTR, GtrModelName)>,
  missing_reason: &str,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let Some(file) = resolved.non_tree_outputs.get(&OutputSelection::Gtr) else {
    return Ok(());
  };
  match model {
    Some((gtr, model_name)) => {
      let gtr_output = GtrOutput::builder().gtr(gtr).model_name(model_name).build();
      json_write_file(&file.path, &gtr_output, JsonPretty(true))
    },
    None => output_unavailable(OutputSelection::Gtr, file, missing_reason, log),
  }
}
