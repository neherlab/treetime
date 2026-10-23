use crate::job::JobId;
use crate::results::outputs::read_auspice;
use crate::results::run_results::finished_record;
use crate::results::tree::{BAD_BRANCH, CATEGORICAL, ResultColoring, result_colorings};
use crate::runs::errors::not_found;
use crate::runs::manager::RunManager;
use deser::Serialize;
use eyre::{Report, WrapErr};
use itertools::izip;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use treetime_io::auspice_types::{AuspiceColoring, AuspiceTree};

const GENOTYPE: &str = "gt";

const MUTED_PALETTE: [&str; 9] = [
  "#332288", "#88ccee", "#44aa99", "#117733", "#999933", "#ddcc77", "#cc6677", "#882255", "#aa4499",
];

pub fn run_auspice(manager: &RunManager, id: &JobId) -> Result<AuspiceDocument, Report> {
  let record = finished_record(manager, id)?;
  let auspice = read_auspice(&manager.store().out_dir(id), &record.output_files)
    .wrap_err_with(|| format!("When reading the Auspice file of run `{}`", id.as_str()))?
    .ok_or_else(|| not_found(format!("run `{}` wrote no Auspice file", id.as_str())))?;
  display_auspice(auspice)
}

pub fn display_auspice(mut auspice: AuspiceTree) -> Result<AuspiceDocument, Report> {
  let colorings = result_colorings(&auspice)?;
  for (coloring, ResultColoring { states, .. }) in izip!(&mut auspice.data.meta.colorings, colorings) {
    if let Some(scale) = muted_scale(coloring, &states) {
      coloring.scale = scale;
    }
  }
  Ok(AuspiceDocument(auspice))
}

/// Auspice JSON of a run, with the color scales the app displays.
#[derive(Clone, Debug, JsonSchema, Serialize)]
#[schemars(transparent)]
pub struct AuspiceDocument(#[schemars(with = "Map<String, Value>")] pub AuspiceTree);

fn muted_scale(coloring: &AuspiceColoring, states: &[String]) -> Option<Vec<[String; 2]>> {
  let muted = coloring.type_ == CATEGORICAL
    && coloring.key != GENOTYPE
    && coloring.key != BAD_BRANCH
    && (1..=MUTED_PALETTE.len()).contains(&states.len());
  muted.then(|| {
    izip!(states, MUTED_PALETTE)
      .map(|(state, color)| [state.clone(), color.to_owned()])
      .collect()
  })
}
