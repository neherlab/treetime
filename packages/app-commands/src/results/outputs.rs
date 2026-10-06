use crate::command::{AppCommand, OutputFile};
use crate::results::homoplasy::HomoplasyStatsFile;
use crate::results::tree::ResultTree;
use app_output::output_plan::OutputSelection;
use app_output::table_output::table_read_file;
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::path::Path;
use treetime::clock::clock_model::ClockModel;
use treetime::clock::rtt::ClockRegressionResult;
use treetime::gtr::get_gtr::GtrOutput;
use treetime::timetree::coalescent::CoalescentSegmentRow;
use treetime::timetree::convergence::metrics::ConvergenceMetrics;
use treetime_io::auspice_types::AuspiceTree;
use treetime_utils::io::json::json_read_file;
use util_augur_node_data_json::AugurNodeDataJsonRefine;

pub fn read_auspice(out_dir: &Path, output_files: &[OutputFile]) -> Result<Option<AuspiceTree>, Report> {
  output_files
    .iter()
    .find(|file| file.kind == OutputSelection::Auspice)
    .map(|file| json_read_file(out_dir.join(&file.path)))
    .transpose()
}

/// An output file of a run that could not be read.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
pub struct OutputProblem {
  /// Path of the file relative to the run's `out/` folder.
  pub path: String,
  /// Why the file could not be read.
  pub message: String,
}

#[derive(Default)]
pub struct RunOutputs {
  pub auspice: Option<AuspiceTree>,
  pub tree: Option<ResultTree>,
  pub clock_model: Option<ClockModel>,
  pub clock_rows: Option<Vec<ClockRegressionResult>>,
  pub node_data: Option<AugurNodeDataJsonRefine>,
  pub gtr: Option<GtrOutput>,
  pub trace: Option<Vec<ConvergenceMetrics>>,
  pub coalescent: Option<Vec<CoalescentSegmentRow>>,
  pub homoplasy: Option<HomoplasyStatsFile>,
  pub problems: Vec<OutputProblem>,
}

impl RunOutputs {
  pub fn read(command: AppCommand, out_dir: &Path, output_files: &[OutputFile]) -> Self {
    let mut outputs = Self::default();
    let wanted = result_outputs(command);
    for file in output_files.iter().filter(|file| wanted.contains(&file.kind)) {
      let path = out_dir.join(&file.path);
      if let Err(err) = outputs.read_file(file.kind, &path) {
        outputs.problems.push(OutputProblem {
          path: file.path.to_string_lossy().into_owned(),
          message: format!("{err:#}"),
        });
      }
    }
    outputs
  }

  fn read_file(&mut self, kind: OutputSelection, path: &Path) -> Result<(), Report> {
    match kind {
      OutputSelection::Auspice => {
        let auspice = json_read_file(path)?;
        self.tree = Some(ResultTree::from_auspice(&auspice)?);
        self.auspice = Some(auspice);
      },
      OutputSelection::ClockModel => self.clock_model = Some(json_read_file(path)?),
      OutputSelection::ClockCsv => self.clock_rows = Some(table_read_file(kind, path)?),
      OutputSelection::AugurNodeData => self.node_data = Some(json_read_file(path)?),
      OutputSelection::Gtr => self.gtr = Some(json_read_file(path)?),
      OutputSelection::Tracelog => self.trace = Some(table_read_file(kind, path)?),
      OutputSelection::CoalescentTsv => self.coalescent = Some(table_read_file(kind, path)?),
      OutputSelection::HomoplasyStats => self.homoplasy = Some(json_read_file(path)?),
      _ => {},
    }
    Ok(())
  }
}

fn result_outputs(command: AppCommand) -> &'static [OutputSelection] {
  match command {
    AppCommand::Timetree => &[
      OutputSelection::Auspice,
      OutputSelection::ClockModel,
      OutputSelection::ClockCsv,
      OutputSelection::AugurNodeData,
      OutputSelection::Tracelog,
      OutputSelection::CoalescentTsv,
    ],
    AppCommand::Clock => &[
      OutputSelection::Auspice,
      OutputSelection::ClockModel,
      OutputSelection::ClockCsv,
    ],
    AppCommand::Optimize => &[
      OutputSelection::Auspice,
      OutputSelection::AugurNodeData,
      OutputSelection::Gtr,
    ],
    AppCommand::Homoplasy => &[OutputSelection::Auspice, OutputSelection::HomoplasyStats],
    AppCommand::Ancestral | AppCommand::Mugration | AppCommand::Prune => &[OutputSelection::Auspice],
  }
}
