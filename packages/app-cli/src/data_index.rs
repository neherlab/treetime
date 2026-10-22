use crate::cli::diagnostics::entry::check_pipeline;
use crate::cli::pipeline::resolve::PipelineDoc;
use app_backend::datasets::{DataIndex, IndexedStep};
use app_commands::command::AppCommand;
use app_commands::config::source::{ConfigSource, parse_config_document};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::fs;
use std::path::Path;
use std::str::FromStr;
use treetime_utils::io::fs::read_file_to_string;
use treetime_utils::make_report;

const PIPELINE_SOURCE: &str = "pipeline.yaml";

pub fn generate_data_index(data_dir: &Path, file_list: &str) -> Result<DataIndex, Report> {
  let paths = file_list
    .split('\0')
    .filter(|path| !path.is_empty())
    .map(str::to_owned)
    .collect_vec();
  DataIndex::new(
    paths,
    |path| file_size(&data_dir.join(path)),
    |path| read_file_to_string(data_dir.join(path)),
    pipeline_steps,
  )
}

pub fn pipeline_steps(text: &str) -> Result<Vec<IndexedStep>, Report> {
  let source = ConfigSource::new(PIPELINE_SOURCE, text);
  let value = parse_config_document(&source, text)?;
  check_pipeline(&source, &value)?;
  PipelineDoc::from_value(value)?
    .steps
    .into_iter()
    .map(|step| {
      let command = AppCommand::from_str(&step.tag)
        .map_err(|err| make_report!("pipeline step `{}` runs an unknown command: {err}", step.name))?;
      Ok(IndexedStep {
        name: step.name,
        command,
      })
    })
    .try_collect()
}

fn file_size(file: &Path) -> Result<u64, Report> {
  fs::metadata(file)
    .map(|metadata| metadata.len())
    .wrap_err_with(|| format!("When reading the size of '{}'", file.display()))
}
