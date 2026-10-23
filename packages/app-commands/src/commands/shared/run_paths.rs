use crate::command::AppCommand;
use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::dependent_settings::setting_flag;
use crate::config::properties::{PathRole, setting_properties};
use app_output::output_plan::{OutputSelection, PathPlaceholder, ResolvedOutputs};
use eyre::Report;
use itertools::Itertools;
use ndarray::Array2;
use serde_json::{Map, Value};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Component, Path, PathBuf, is_separator};
use treetime_utils::io::file::{is_path_stdin, is_path_stdout};
use treetime_utils::io::fs::absolute_path;

const PLACEHOLDERS: [PathPlaceholder; 2] = [PathPlaceholder::Alignment, PathPlaceholder::Cds];

const OUTPUT_ALL_KEY: &str = "output_all";

pub fn check_command_paths(
  command: AppCommand,
  settings: &Map<String, Value>,
  outputs: &ResolvedOutputs,
) -> Result<(), Report> {
  let inputs = run_inputs(command, settings)?;
  let outputs = run_outputs(settings, &outputs.paths_by_selection());
  check_run_paths(&inputs, &outputs)
}

pub fn check_run_paths(inputs: &[RunPath], outputs: &[RunPath]) -> Result<(), Report> {
  reject_repeated_standard_input(inputs)?;
  let outputs = outputs
    .iter()
    .filter(|output| !output.is_stdio())
    .map(|output| Ok((output, output.pattern()?)))
    .collect::<Result<Vec<_>, Report>>()?;
  for input in inputs.iter().filter(|input| !input.is_stdio()) {
    let input_pattern = input.pattern()?;
    if let Some((output, output_pattern)) = outputs
      .iter()
      .find(|(_, output_pattern)| patterns_intersect(&input_pattern, output_pattern))
    {
      let verb = if is_literal(&input_pattern) && is_literal(output_pattern) {
        "names"
      } else {
        "can name"
      };
      return Err(argument_error(
        &[input.setting.as_str(), output.setting.as_str()],
        format!(
          "{} '{}' {verb} the same file as the output {} '{}': the run would overwrite its own input",
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

pub fn run_outputs(
  settings: &Map<String, Value>,
  by_selection: &BTreeMap<OutputSelection, Vec<PathBuf>>,
) -> Vec<RunPath> {
  by_selection
    .iter()
    .flat_map(|(selection, paths)| {
      let key = selection.setting_key();
      let setting = if settings.get(&key).is_some_and(|value| !value.is_null()) {
        key
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

  fn pattern(&self) -> Result<Vec<PatternChar>, Report> {
    let path = opened_path(&absolute_path(&self.path)?);
    let text = path.to_string_lossy();
    Ok(match self.role {
      PathRole::Input => text.chars().map(PatternChar::Literal).collect(),
      PathRole::InputTemplate | PathRole::Output => template_pattern(&text),
    })
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PatternChar {
  Literal(char),
  AnyInName,
}

impl PatternChar {
  fn is_separator(self) -> bool {
    matches!(self, Self::Literal(c) if is_separator(c))
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

fn setting_paths(value: Option<&Value>) -> Vec<String> {
  match value {
    Some(Value::String(path)) => vec![path.clone()],
    Some(Value::Array(items)) => items.iter().filter_map(Value::as_str).map(ToOwned::to_owned).collect(),
    _ => vec![],
  }
}

fn template_pattern(text: &str) -> Vec<PatternChar> {
  let tokens = PLACEHOLDERS
    .iter()
    .flat_map(|placeholder| placeholder.tokens())
    .collect_vec();
  let mut pattern = Vec::with_capacity(text.len());
  let mut rest = text;
  while !rest.is_empty() {
    if let Some(after) = tokens.iter().find_map(|token| rest.strip_prefix(*token)) {
      pattern.push(PatternChar::AnyInName);
      rest = after;
    } else {
      let mut chars = rest.chars();
      pattern.extend(chars.next().map(PatternChar::Literal));
      rest = chars.as_str();
    }
  }
  pattern
}

fn is_literal(pattern: &[PatternChar]) -> bool {
  !pattern.contains(&PatternChar::AnyInName)
}

fn patterns_intersect(left: &[PatternChar], right: &[PatternChar]) -> bool {
  let (left_len, right_len) = (left.len(), right.len());
  let mut matches = Array2::from_elem((left_len + 1, right_len + 1), false);
  matches[[left_len, right_len]] = true;
  for i in (0..=left_len).rev() {
    for j in (0..=right_len).rev() {
      if i == left_len && j == right_len {
        continue;
      }
      matches[[i, j]] = match (left.get(i), right.get(j)) {
        (Some(PatternChar::AnyInName), next) => {
          matches[[i + 1, j]] || (next.is_some_and(|c| !c.is_separator()) && matches[[i, j + 1]])
        },
        (next, Some(PatternChar::AnyInName)) => {
          matches[[i, j + 1]] || (next.is_some_and(|c| !c.is_separator()) && matches[[i + 1, j]])
        },
        (Some(PatternChar::Literal(left_char)), Some(PatternChar::Literal(right_char))) => {
          left_char == right_char && matches[[i + 1, j + 1]]
        },
        _ => false,
      };
    }
  }
  matches[[0, 0]]
}

fn and_list(items: &[String]) -> String {
  match items.split_last() {
    Some((last, init)) if !init.is_empty() => format!("{} and {last}", init.join(", ")),
    Some((last, _)) => last.clone(),
    None => String::new(),
  }
}
