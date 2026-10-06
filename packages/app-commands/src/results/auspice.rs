use crate::job::JobId;
use crate::results::outputs::read_auspice;
use crate::results::run_results::finished_record;
use crate::results::tree::{ResultColoring, result_colorings};
use crate::runs::errors::not_found;
use crate::runs::manager::RunManager;
use eyre::{Report, WrapErr};
use itertools::izip;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use treetime_io::auspice_types::AuspiceTree;

/// Auspice JSON of a run, with the color scales the app displays.
#[derive(Clone, Debug, JsonSchema, deser::Serialize, deser::Deserialize)]
#[schemars(transparent)]
pub struct AuspiceDocument(#[schemars(with = "Map<String, Value>")] pub AuspiceTree);

pub fn run_auspice(manager: &RunManager, id: &JobId) -> Result<AuspiceDocument, Report> {
  let record = finished_record(manager, id)?;
  let auspice = read_auspice(&manager.store().out_dir(id), &record.output_files)
    .wrap_err_with(|| format!("When reading the Auspice file of run `{}`", id.as_str()))?
    .ok_or_else(|| not_found(format!("run `{}` wrote no Auspice file", id.as_str())))?;
  display_auspice(auspice)
}

pub fn display_auspice(mut auspice: AuspiceTree) -> Result<AuspiceDocument, Report> {
  let colorings = result_colorings(&auspice)?;
  for (coloring, ResultColoring { scale, .. }) in izip!(&mut auspice.data.meta.colorings, colorings) {
    coloring.scale = scale.into_iter().map(|color| [color.state, color.color]).collect();
  }
  Ok(AuspiceDocument(auspice))
}
