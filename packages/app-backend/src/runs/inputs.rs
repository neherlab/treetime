use crate::runs::input_store::InputStore;
use crate::runs::record::RunInput;
use crate::runs::setting_differences::ComparedInput;
use app_commands::command::AppCommand;
use app_commands::commands::ancestral::aa_node_data::{cds_names, translation_path};
use app_commands::config::properties::{PathRole, setting_properties};
use deser::adapters::As;
use deser_serde::Serde;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::warn;
use parking_lot::Mutex;
use serde_json::{Map, Value};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
use std::fmt::{self, Write};
use std::fs::{self, File, Metadata};
use std::io::{self, Read};
#[cfg(unix)]
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};
use std::time::SystemTime;
use treetime_utils::io::json::{JsonPretty, json_write_str};

const HASH_BUFFER_SIZE: usize = 1 << 16;

pub fn hash_inputs(command: AppCommand, config: &Map<String, Value>) -> Result<HashedInputs, Report> {
  let mut canonical = config.clone();
  let mut inputs = vec![];
  for property in setting_properties(command.config_schema().as_value())? {
    let setting = property.key.as_str();
    match property.path_role {
      Some(PathRole::Output) => {
        canonical.remove(setting);
      },
      Some(PathRole::Input) => {
        if let Some(value) = canonical.get_mut(setting) {
          hash_path_values(setting, value, &mut inputs)?;
        }
      },
      Some(PathRole::InputTemplate) => {
        let paths = template_input_paths(config, setting)?;
        if let Some(value) = canonical.get_mut(setting) {
          let recorded: Vec<RunInput> = paths.iter().map(|path| record_input(setting, path)).try_collect()?;
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
  let canonical = json_write_str(
    &As::<_, Serde>::new(&sorted_keys(&Value::Object(canonical))),
    JsonPretty(false),
  )?;
  Ok(HashedInputs {
    inputs,
    config_hash: sha256_hex(canonical.as_bytes())?,
  })
}

pub fn draft_inputs(
  command: AppCommand,
  config: &Map<String, Value>,
  sha256: impl Fn(&Path) -> Result<String, Report>,
) -> Result<Vec<ComparedInput>, Report> {
  let mut inputs = vec![];
  for property in setting_properties(command.config_schema().as_value())? {
    let setting = property.key.as_str();
    let paths = match property.path_role {
      Some(PathRole::Input) => Ok(config.get(setting).map_or_else(Vec::new, path_values)),
      Some(PathRole::InputTemplate) => template_input_paths(config, setting),
      Some(PathRole::Output) | None => continue,
    };
    match paths {
      Ok(paths) => inputs.extend(paths.into_iter().map(|path| {
        let sha256 = sha256(&path)
          .inspect_err(|report| {
            warn!(
              "The input '{}' of setting `{setting}` counts as changed: {report:#}",
              path.display()
            );
          })
          .ok();
        ComparedInput {
          setting: setting.to_owned(),
          path,
          sha256,
        }
      })),
      Err(report) => {
        warn!("The inputs of setting `{setting}` count as changed: {report:#}");
        inputs.push(ComparedInput {
          setting: setting.to_owned(),
          path: config
            .get(setting)
            .and_then(Value::as_str)
            .map(PathBuf::from)
            .unwrap_or_default(),
          sha256: None,
        });
      },
    }
  }
  Ok(inputs)
}

/// SHA-256 of input files, remembered for files whose size, modification time, and identity stay the same.
#[derive(Default)]
pub struct InputHashes {
  known: Mutex<BTreeMap<PathBuf, KnownHash>>,
}

impl InputHashes {
  pub fn sha256(&self, store: &InputStore, path: &Path) -> Result<String, Report> {
    if let Some(key) = store.key_of(path) {
      return Ok(key.sha256);
    }
    let metadata = fs::metadata(path).wrap_err_with(|| format!("When reading the size of '{}'", path.display()))?;
    let Ok(stamp) = FileStamp::of(&metadata) else {
      return Ok(file_sha256(path)?.1);
    };
    if let Some(known) = self.known.lock().get(path).filter(|known| known.stamp == stamp) {
      return Ok(known.sha256.clone());
    }
    let (_, sha256) = file_sha256(path)?;
    self.known.lock().insert(
      path.to_path_buf(),
      KnownHash {
        stamp,
        sha256: sha256.clone(),
      },
    );
    Ok(sha256)
  }
}

struct KnownHash {
  stamp: FileStamp,
  sha256: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct FileStamp {
  size: u64,
  modified: SystemTime,
  identity: FileIdentity,
}

impl FileStamp {
  fn of(metadata: &Metadata) -> io::Result<Self> {
    Ok(Self {
      size: metadata.len(),
      modified: metadata.modified()?,
      identity: FileIdentity::of(metadata)?,
    })
  }
}

#[cfg(unix)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct FileIdentity {
  device: u64,
  inode: u64,
  changed_seconds: i64,
  changed_nanoseconds: i64,
}

#[cfg(unix)]
impl FileIdentity {
  fn of(metadata: &Metadata) -> io::Result<Self> {
    Ok(Self {
      device: metadata.dev(),
      inode: metadata.ino(),
      changed_seconds: metadata.ctime(),
      changed_nanoseconds: metadata.ctime_nsec(),
    })
  }
}

#[cfg(not(unix))]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct FileIdentity {
  created: SystemTime,
}

#[cfg(not(unix))]
impl FileIdentity {
  fn of(metadata: &Metadata) -> io::Result<Self> {
    Ok(Self {
      created: metadata.created()?,
    })
  }
}

pub struct HashedInputs {
  pub inputs: Vec<RunInput>,
  pub config_hash: String,
}

pub fn sha256_hex(bytes: &[u8]) -> Result<String, Report> {
  Ok(hex_digest(&Sha256::digest(bytes))?)
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

fn path_values(value: &Value) -> Vec<PathBuf> {
  match value {
    Value::String(path) => vec![PathBuf::from(path)],
    Value::Array(items) => items.iter().flat_map(path_values).collect(),
    _ => vec![],
  }
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
  Ok((size, hex_digest(&hasher.finalize())?))
}

pub fn hex_digest(bytes: &[u8]) -> Result<String, fmt::Error> {
  bytes.iter().try_fold(String::new(), |mut hex, byte| {
    write!(hex, "{byte:02x}")?;
    Ok(hex)
  })
}

pub fn template_input_paths(config: &Map<String, Value>, key: &str) -> Result<Vec<PathBuf>, Report> {
  let Some(Value::String(template)) = config.get(key) else {
    return Ok(vec![]);
  };
  Ok(
    cds_names(config)?
      .iter()
      .map(|cds| translation_path(template, cds))
      .collect(),
  )
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
