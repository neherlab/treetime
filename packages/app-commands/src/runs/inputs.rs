use crate::command::AppCommand;
use crate::commands::ancestral::aa_node_data::translation_input_paths;
use crate::config::properties::{PathRole, leaf_properties};
use crate::config::settings::{remove_setting, setting_mut, setting_ref};
use crate::runs::record::RunInput;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use serde_json::{Map, Value};
use sha2::{Digest, Sha256};
use std::fmt::{self, Write};
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

const HASH_BUFFER_SIZE: usize = 1 << 16;

pub struct HashedInputs {
  pub inputs: Vec<RunInput>,
  pub config_hash: String,
}

pub fn hash_inputs(command: AppCommand, config: &Map<String, Value>) -> Result<HashedInputs, Report> {
  let mut canonical = config.clone();
  let mut inputs = vec![];
  for leaf in leaf_properties(command.config_schema().as_value())? {
    let setting = leaf.key_path.join(".");
    match leaf.path_role {
      Some(PathRole::Output) => remove_setting(&mut canonical, &leaf.key_path),
      Some(PathRole::Input) => {
        if let Some(value) = setting_mut(&mut canonical, &leaf.key_path) {
          hash_path_values(&setting, value, &mut inputs)?;
        }
      },
      Some(PathRole::InputTemplate) => {
        let paths = template_inputs(config, &leaf.key_path)?;
        if let Some(value) = setting_mut(&mut canonical, &leaf.key_path) {
          let recorded: Vec<RunInput> = paths.iter().map(|path| record_input(&setting, path)).try_collect()?;
          *value = Value::Array(
            recorded
              .iter()
              .map(|input| Value::String(input.sha256.clone()))
              .collect(),
          );
          inputs.extend(recorded);
        }
      },
      None => {},
    }
  }
  let canonical = serde_json::to_string(&sorted_keys(&Value::Object(canonical)))?;
  Ok(HashedInputs {
    inputs,
    config_hash: sha256_hex(canonical.as_bytes())?,
  })
}

pub fn file_sha256(path: &Path) -> Result<(usize, String), Report> {
  let mut file = File::open(path).wrap_err_with(|| format!("When opening input '{}'", path.display()))?;
  let mut hasher = Sha256::new();
  let mut buffer = vec![0_u8; HASH_BUFFER_SIZE];
  let mut size: usize = 0;
  loop {
    let read = file
      .read(&mut buffer)
      .wrap_err_with(|| format!("When reading input '{}'", path.display()))?;
    if read == 0 {
      break;
    }
    hasher.update(&buffer[..read]);
    size += read;
  }
  Ok((size, hex(&hasher.finalize())?))
}

pub fn sha256_hex(bytes: &[u8]) -> Result<String, Report> {
  Ok(hex(&Sha256::digest(bytes))?)
}

fn hex(bytes: &[u8]) -> Result<String, fmt::Error> {
  bytes.iter().try_fold(String::new(), |mut hex, byte| {
    write!(hex, "{byte:02x}")?;
    Ok(hex)
  })
}

fn hash_path_values(setting: &str, value: &mut Value, inputs: &mut Vec<RunInput>) -> Result<(), Report> {
  match value {
    Value::String(path) => {
      let input = record_input(setting, Path::new(path))?;
      *value = Value::String(input.sha256.clone());
      inputs.push(input);
    },
    Value::Array(items) => {
      for item in items {
        hash_path_values(setting, item, inputs)?;
      }
    },
    _ => {},
  }
  Ok(())
}

fn record_input(setting: &str, path: &Path) -> Result<RunInput, Report> {
  let (size, sha256) = file_sha256(path).wrap_err_with(|| format!("When hashing the input of setting `{setting}`"))?;
  Ok(RunInput {
    setting: setting.to_owned(),
    path: path.to_path_buf(),
    size,
    sha256,
  })
}

fn template_inputs(config: &Map<String, Value>, key_path: &[String]) -> Result<Vec<PathBuf>, Report> {
  let Some(Value::String(template)) = setting_ref(config, key_path) else {
    return Ok(vec![]);
  };
  let cdses = config
    .get("cdses")
    .and_then(Value::as_array)
    .map(|names| names.iter().filter_map(Value::as_str).map(str::to_owned).collect_vec())
    .unwrap_or_default();
  let annotation = config.get("annotation").and_then(Value::as_str).map(Path::new);
  translation_input_paths(template, &cdses, annotation)
}

fn sorted_keys(value: &Value) -> Value {
  match value {
    Value::Object(object) => Value::Object(
      object
        .iter()
        .sorted_by(|(a, _), (b, _)| a.cmp(b))
        .map(|(key, value)| (key.clone(), sorted_keys(value)))
        .collect::<Map<String, Value>>(),
    ),
    Value::Array(items) => Value::Array(items.iter().map(sorted_keys).collect()),
    other => other.clone(),
  }
}
