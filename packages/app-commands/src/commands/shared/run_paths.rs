use crate::command::AppCommand;
use crate::commands::ancestral::aa_node_data::{cds_output_paths, translation_path};
use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::dependent_settings::setting_flag;
use crate::commands::shared::output_args::OutputRequest;
use crate::config::properties::{PathRole, setting_properties};
use app_output::output_plan::{OutputSelection, ResolvedOutputs};
use eyre::Report;
use itertools::Itertools;
use serde_json::{Map, Value};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Component, Path, PathBuf};
use treetime_utils::io::file::{is_path_stdin, is_path_stdout};
use treetime_utils::io::fs::absolute_path;

const OUTPUT_ALL_KEY: &str = "output_all";

const CDS_SELECTION_KEY: &str = "cds_selection";

pub(crate) fn check_command_paths(
  command: AppCommand,
  settings: &Map<String, Value>,
  request: &OutputRequest,
) -> Result<(), Report> {
  let paths = RunPaths::new(run_inputs(command, settings)?, request, &request.resolve()?);
  paths.check_without_cds_names()?;
  let cdses = listed_cds_names(settings);
  if cdses.is_empty() {
    return Ok(());
  }
  paths.check_with_cds_names(&cdses)
}

pub fn run_inputs(command: AppCommand, settings: &Map<String, Value>) -> Result<Vec<RunPath>, Report> {
  let schema = command.config_schema();
  Ok(
    setting_properties(schema.as_value())?
      .into_iter()
      .filter_map(|property| match property.path_role? {
        role @ (PathRole::Input | PathRole::InputTemplate) => Some((property.key, role)),
        PathRole::Output => None,
      })
      .flat_map(|(setting, role)| {
        setting_paths(settings.get(&setting))
          .into_iter()
          .map(move |path| RunPath {
            setting: setting.clone(),
            path,
            role,
          })
      })
      .collect(),
  )
}

pub(crate) fn run_outputs(
  request: &OutputRequest,
  by_selection: &BTreeMap<OutputSelection, Vec<PathBuf>>,
) -> Vec<RunPath> {
  let given = request.given_selections();
  by_selection
    .iter()
    .flat_map(|(selection, paths)| {
      let setting = if given.contains(selection) {
        selection.setting_key()
      } else {
        OUTPUT_ALL_KEY.to_owned()
      };
      paths.iter().map(move |path| RunPath {
        setting: setting.clone(),
        path: path.to_string_lossy().into_owned(),
        role: PathRole::Output,
      })
    })
    .collect()
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RunPaths {
  pub(crate) inputs: Vec<RunPath>,
  pub(crate) outputs: Vec<RunPath>,
  pub(crate) cds_outputs: Vec<RunPath>,
}

impl RunPaths {
  pub(crate) fn new(inputs: Vec<RunPath>, request: &OutputRequest, resolved: &ResolvedOutputs) -> Self {
    let (cds_outputs, outputs): (BTreeMap<_, _>, BTreeMap<_, _>) = resolved
      .paths_by_selection()
      .into_iter()
      .partition(|(selection, _)| *selection == OutputSelection::Translations);
    Self {
      inputs,
      outputs: run_outputs(request, &outputs),
      cds_outputs: run_outputs(request, &cds_outputs),
    }
  }

  pub(crate) fn check_without_cds_names(&self) -> Result<(), Report> {
    reject_repeated_standard_input(&self.inputs)?;
    reject_overwritten_inputs(
      self.inputs.iter().filter(|input| input.role != PathRole::InputTemplate),
      &self.outputs,
    )
  }

  pub(crate) fn check_with_cds_names(&self, cdses: &[String]) -> Result<(), Report> {
    let inputs = self
      .inputs
      .iter()
      .flat_map(|input| {
        if input.role == PathRole::InputTemplate {
          cdses
            .iter()
            .map(|cds| input.expanded(&translation_path(&input.path, cds)))
            .collect_vec()
        } else {
          vec![input.clone()]
        }
      })
      .collect_vec();
    let cds_outputs = self
      .cds_outputs
      .iter()
      .map(|output| {
        Ok(
          cds_output_paths(Path::new(&output.path), cdses)?
            .into_values()
            .map(|path| output.expanded(&path))
            .collect_vec(),
        )
      })
      .flatten_ok()
      .collect::<Result<Vec<_>, Report>>()?;
    reject_overwritten_inputs(&inputs, self.outputs.iter().chain(&cds_outputs))
  }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RunPath {
  pub setting: String,
  pub path: String,
  pub role: PathRole,
}

impl RunPath {
  pub fn is_stdio(&self) -> bool {
    match self.role {
      PathRole::Input | PathRole::InputTemplate => is_path_stdin(&self.path),
      PathRole::Output => is_path_stdout(&self.path),
    }
  }

  fn expanded(&self, path: &Path) -> Self {
    Self {
      setting: self.setting.clone(),
      path: path.to_string_lossy().into_owned(),
      role: match self.role {
        PathRole::Input | PathRole::InputTemplate => PathRole::Input,
        PathRole::Output => PathRole::Output,
      },
    }
  }

  fn opened(&self) -> Result<PathBuf, Report> {
    Ok(opened_path(&absolute_path(&self.path)?))
  }
}

fn reject_repeated_standard_input(inputs: &[RunPath]) -> Result<(), Report> {
  let readers = inputs
    .iter()
    .filter(|input| input.is_stdio())
    .map(|input| input.setting.as_str())
    .collect_vec();
  if readers.len() < 2 {
    return Ok(());
  }
  let settings = readers.iter().copied().unique().collect_vec();
  let message = match settings.as_slice() {
    [setting] => format!(
      "Standard input can be read only once, but {} lists it {} times",
      setting_flag(setting),
      readers.len()
    ),
    _ => format!(
      "Standard input can be read only once, but {} read it",
      and_list(&settings.iter().map(|setting| setting_flag(setting)).collect_vec())
    ),
  };
  Err(argument_error(&settings, message))
}

fn reject_overwritten_inputs<'a>(
  inputs: impl IntoIterator<Item = &'a RunPath>,
  outputs: impl IntoIterator<Item = &'a RunPath>,
) -> Result<(), Report> {
  let outputs = outputs
    .into_iter()
    .filter(|output| !output.is_stdio())
    .map(|output| Ok((output, output.opened()?)))
    .collect::<Result<Vec<_>, Report>>()?;
  for input in inputs.into_iter().filter(|input| !input.is_stdio()) {
    let opened = input.opened()?;
    if let Some((output, _)) = outputs.iter().find(|(_, output_opened)| *output_opened == opened) {
      return Err(argument_error(
        &[input.setting.as_str(), output.setting.as_str()],
        format!(
          "{} '{}' names the same file as the output {} '{}': the run would overwrite its own input",
          setting_flag(&input.setting),
          input.path,
          setting_flag(&output.setting),
          output.path
        ),
      ));
    }
  }
  Ok(())
}

#[cfg_attr(
  dylint_lib = "custom",
  expect(
    error_dropped_by_pattern,
    reason = "a prefix that does not resolve only ends the search at a shorter prefix, and the rest of the path is compared as written"
  )
)]
fn opened_path(path: &Path) -> PathBuf {
  let components = path.components().collect_vec();
  let (existing, rest) = (0..=components.len())
    .rev()
    .find_map(
      |count| match fs::canonicalize(components[..count].iter().collect::<PathBuf>()) {
        Ok(canonical) => Some((canonical, &components[count..])),
        Err(_) => None,
      },
    )
    .unwrap_or((PathBuf::new(), &components[..]));
  rest.iter().fold(existing, |mut path, component| {
    match component {
      Component::ParentDir => {
        path.pop();
      },
      Component::CurDir => {},
      Component::Prefix(_) | Component::RootDir | Component::Normal(_) => path.push(component),
    }
    path
  })
}

fn listed_cds_names(settings: &Map<String, Value>) -> Vec<String> {
  settings
    .get(CDS_SELECTION_KEY)
    .and_then(Value::as_array)
    .map_or_else(Vec::new, |names| {
      names.iter().filter_map(Value::as_str).map(ToOwned::to_owned).collect()
    })
}

fn setting_paths(value: Option<&Value>) -> Vec<String> {
  match value {
    Some(Value::String(path)) => vec![path.clone()],
    Some(Value::Array(items)) => items.iter().filter_map(Value::as_str).map(ToOwned::to_owned).collect(),
    _ => vec![],
  }
}

fn and_list(items: &[String]) -> String {
  match items.split_last() {
    Some((last, init)) if !init.is_empty() => format!("{} and {last}", init.join(", ")),
    Some((last, _)) => last.clone(),
    None => String::new(),
  }
}
