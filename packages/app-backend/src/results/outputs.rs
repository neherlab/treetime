use crate::command::OutputFile;
use crate::results::homoplasy::HomoplasyStatsFile;
use crate::results::tree::ResultTree;
use app_commands::command::AppCommand;
use app_output::output_plan::OutputSelection;
use app_output::table_output::table_read_file;
use deser::Serialize;
use eyre::Report;
use schemars::JsonSchema;
use std::path::Path;
use treetime::clock::clock_model::ClockModel;
use treetime::clock::rtt::ClockRegressionResult;
use treetime::gtr::get_gtr::ModelOutput;
use treetime::timetree::coalescent::CoalescentSegmentRow;
use treetime::timetree::convergence::metrics::ConvergenceMetrics;
use treetime_io::auspice_types::AuspiceTree;
use treetime_utils::io::json::json_read_file;
use util_augur_node_data_json::AugurNodeDataJsonRefine;

pub fn read_auspice(out_dir: &Path, output_files: &[OutputFile]) -> Result<Option<AuspiceTree>, Report> {
  shown_auspice(output_files)
    .map(|file| json_read_file(out_dir.join(&file.path)))
    .transpose()
}

#[derive(Default)]
pub struct RunOutputs {
  pub auspice: Option<AuspiceTree>,
  pub tree: Option<ResultTree>,
  pub clock_model: Option<ClockModel>,
  pub clock_rows: Option<Vec<ClockRegressionResult>>,
  pub node_data: Option<AugurNodeDataJsonRefine>,
  pub gtr: Option<ModelOutput>,
  pub trace: Option<Vec<ConvergenceMetrics>>,
  pub coalescent: Option<Vec<CoalescentSegmentRow>>,
  pub homoplasy: Option<HomoplasyStatsFile>,
  pub problems: Vec<OutputProblem>,
}

impl RunOutputs {
  pub fn read(command: AppCommand, out_dir: &Path, output_files: &[OutputFile]) -> Self {
    let mut outputs = Self::default();
    let wanted = result_outputs(command);
    let auspice = shown_auspice(output_files);
    let read = output_files
      .iter()
      .filter(|file| wanted.contains(&file.kind) && (file.kind != OutputSelection::Auspice || Some(*file) == auspice));
    for file in read {
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
      OutputSelection::Rtt => self.clock_rows = Some(table_read_file(kind, path)?),
      OutputSelection::NodeData => self.node_data = Some(json_read_file(path)?),
      OutputSelection::Model => self.gtr = Some(json_read_file(path)?),
      OutputSelection::Tracelog => self.trace = Some(table_read_file(kind, path)?),
      OutputSelection::Coalescent => self.coalescent = Some(table_read_file(kind, path)?),
      OutputSelection::HomoplasyStats => self.homoplasy = Some(json_read_file(path)?),
      _ => {},
    }
    Ok(())
  }
}

pub fn shown_auspice(output_files: &[OutputFile]) -> Option<&OutputFile> {
  output_files
    .iter()
    .filter(|file| file.kind == OutputSelection::Auspice)
    .min_by_key(|file| file.partition.as_ref().map_or(0, |partition| partition.number))
}

/// An output file of a run that could not be read.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct OutputProblem {
  /// Path of the file relative to the run's `out/` folder.
  pub path: String,
  /// Why the file could not be read.
  pub message: String,
}

fn result_outputs(command: AppCommand) -> &'static [OutputSelection] {
  match command {
    AppCommand::Timetree => &[
      OutputSelection::Auspice,
      OutputSelection::ClockModel,
      OutputSelection::Rtt,
      OutputSelection::NodeData,
      OutputSelection::Tracelog,
      OutputSelection::Coalescent,
      OutputSelection::Model,
    ],
    AppCommand::Clock => &[
      OutputSelection::Auspice,
      OutputSelection::ClockModel,
      OutputSelection::Rtt,
    ],
    AppCommand::Optimize => &[
      OutputSelection::Auspice,
      OutputSelection::NodeData,
      OutputSelection::Model,
    ],
    AppCommand::Homoplasy => &[OutputSelection::Auspice, OutputSelection::HomoplasyStats],
    AppCommand::Ancestral | AppCommand::Mugration | AppCommand::Prune => {
      &[OutputSelection::Auspice, OutputSelection::Model]
    },
  }
}
