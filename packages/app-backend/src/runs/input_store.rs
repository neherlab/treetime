use crate::runs::errors::{UploadTooLarge, invalid};
use crate::runs::files::{files_size, relative_name};
use crate::runs::inputs::hex_digest;
use deser::{Deserialize, Serialize};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::JsonSchema;
use sha2::{Digest, Sha256};
use std::collections::BTreeSet;
use std::fs;
use std::io::{self, Read, Write};
use std::path::{Component, Path, PathBuf};
use tempfile::NamedTempFile;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::fmt::float::float_to_significant_digits;

const STORE_DIR: &str = ".inputs";

const SETS_DIR: &str = "sets";

const UPLOAD_BUFFER_SIZE: usize = 1 << 16;

const SHA256_HEX_LEN: usize = 64;

#[derive(Clone, Debug)]
pub struct InputStore {
  root: PathBuf,
}

impl InputStore {
  pub fn new(runs_root: &Path) -> Self {
    Self {
      root: runs_root.join(STORE_DIR),
    }
  }

  pub fn root(&self) -> &Path {
    &self.root
  }

  pub fn path(&self, key: &StoredInput) -> PathBuf {
    self.root.join(&key.sha256).join(&key.name)
  }

  pub fn key_of(&self, path: &Path) -> Option<StoredInput> {
    let mut components = path.strip_prefix(&self.root).ok()?.components();
    let (Some(Component::Normal(sha256)), Some(Component::Normal(name)), None) =
      (components.next(), components.next(), components.next())
    else {
      return None;
    };
    let sha256 = sha256.to_str().filter(|sha256| is_sha256_hex(sha256))?;
    let name = name.to_str().filter(|name| validate_upload_name(name).is_ok())?;
    Some(StoredInput {
      sha256: sha256.to_owned(),
      name: name.to_owned(),
    })
  }

  pub fn size(&self, key: &StoredInput) -> Result<usize, Report> {
    let path = self.path(key);
    let metadata = fs::metadata(&path).wrap_err_with(|| format!("When reading the size of '{}'", path.display()))?;
    Ok(usize::try_from(metadata.len())?)
  }

  pub fn store(&self, name: &str, body: &mut dyn Read, allowed: usize, limit: usize) -> Result<StoredUpload, Report> {
    let mut writer = self.begin(name, allowed, limit)?;
    let mut buffer = vec![0_u8; UPLOAD_BUFFER_SIZE];
    loop {
      let read = body.read(&mut buffer).wrap_err("When receiving the upload")?;
      if read == 0 {
        return writer.finish();
      }
      writer.write(&buffer[..read])?;
    }
  }

  pub fn begin(&self, name: &str, allowed: usize, limit: usize) -> Result<UploadWriter, Report> {
    validate_upload_name(name)?;
    fs::create_dir_all(&self.root)
      .wrap_err_with(|| format!("When creating the input store '{}'", self.root.display()))?;
    let file = NamedTempFile::new_in(&self.root)
      .wrap_err_with(|| format!("When creating a file in '{}'", self.root.display()))?;
    Ok(UploadWriter {
      file,
      hasher: Sha256::new(),
      name: name.to_owned(),
      size: 0,
      allowed,
      limit,
    })
  }

  pub fn commit(&self, upload: StoredUpload) -> Result<StoredInput, Report> {
    let StoredUpload { file, key, .. } = upload;
    let target = self.path(&key);
    if target.is_file() {
      return Ok(key);
    }
    let folder = self.root.join(&key.sha256);
    fs::create_dir_all(&folder).wrap_err_with(|| format!("When creating '{}'", folder.display()))?;
    file
      .persist(&target)
      .map_err(|err| Report::new(err.error))
      .wrap_err_with(|| format!("When saving the upload '{}'", target.display()))?;
    Ok(key)
  }

  pub fn keys(&self) -> Result<BTreeSet<StoredInput>, Report> {
    let mut keys = BTreeSet::new();
    for folder in read_dir_if_exists(&self.root)? {
      let Some(sha256) = folder.file_name().and_then(|name| name.to_str()).map(str::to_owned) else {
        continue;
      };
      if !is_sha256_hex(&sha256) || !folder.is_dir() {
        continue;
      }
      for file in read_dir_if_exists(&folder)? {
        if let Some(key) = self.key_of(&file) {
          keys.insert(key);
        }
      }
    }
    Ok(keys)
  }

  pub fn remove_partial_uploads(&self) -> Result<(), Report> {
    for path in read_dir_if_exists(&self.root)? {
      if path.is_file() {
        fs::remove_file(&path).wrap_err_with(|| format!("When removing the partial upload '{}'", path.display()))?;
      }
    }
    Ok(())
  }

  pub fn missing(&self, paths: &[PathBuf]) -> Vec<PathBuf> {
    paths
      .iter()
      .filter(|path| self.key_of(path).is_some_and(|key| !self.path(&key).is_file()))
      .cloned()
      .collect()
  }

  pub fn remove(&self, key: &StoredInput) -> Result<(), Report> {
    let path = self.path(key);
    match fs::remove_file(&path) {
      Err(err) if err.kind() != io::ErrorKind::NotFound => {
        return Err(Report::new(err).wrap_err(format!("When removing the stored input '{}'", path.display())));
      },
      _ => {},
    }
    let folder = self.root.join(&key.sha256);
    match fs::remove_dir(&folder) {
      Err(err) if !matches!(err.kind(), io::ErrorKind::NotFound | io::ErrorKind::DirectoryNotEmpty) => {
        Err(Report::new(err).wrap_err(format!("When removing the folder '{}'", folder.display())))
      },
      _ => Ok(()),
    }
  }
  pub fn set_dir(&self, sha256: &str) -> PathBuf {
    self.root.join(SETS_DIR).join(sha256)
  }

  pub fn set_of(&self, path: &Path) -> Option<(String, PathBuf)> {
    let mut components = path.strip_prefix(self.root.join(SETS_DIR)).ok()?.components();
    let Some(Component::Normal(sha256)) = components.next() else {
      return None;
    };
    let sha256 = sha256.to_str().filter(|sha256| is_sha256_hex(sha256))?;
    let rest = components.as_path();
    rest
      .components()
      .all(|component| matches!(component, Component::Normal(_)))
      .then(|| (sha256.to_owned(), rest.to_path_buf()))
  }

  pub fn commit_set(&self, staged: &Path, sha256: &str) -> Result<(), Report> {
    let target = self.set_dir(sha256);
    if target.is_dir() {
      return Ok(());
    }
    let sets = self.root.join(SETS_DIR);
    fs::create_dir_all(&sets).wrap_err_with(|| format!("When creating '{}'", sets.display()))?;
    fs::rename(staged, &target).wrap_err_with(|| format!("When saving the file set '{}'", target.display()))
  }

  pub fn set_keys(&self) -> Result<BTreeSet<String>, Report> {
    Ok(
      read_dir_if_exists(&self.root.join(SETS_DIR))?
        .into_iter()
        .filter(|folder| folder.is_dir())
        .filter_map(|folder| folder.file_name().and_then(|name| name.to_str()).map(str::to_owned))
        .filter(|sha256| is_sha256_hex(sha256))
        .collect(),
    )
  }

  pub fn set_size(&self, sha256: &str) -> Result<usize, Report> {
    files_size(&self.set_dir(sha256))
  }

  pub fn remove_set(&self, sha256: &str) -> Result<(), Report> {
    let folder = self.set_dir(sha256);
    match fs::remove_dir_all(&folder) {
      Err(err) if err.kind() != io::ErrorKind::NotFound => {
        Err(Report::new(err).wrap_err(format!("When removing the file set '{}'", folder.display())))
      },
      _ => Ok(()),
    }
  }
}

/// Input paths, for example those of a form, to look up in the store of uploaded inputs.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct StoredInputsRequest {
  /// Paths of the input files, as the form holds them.
  pub paths: Vec<PathBuf>,
}

/// The inputs of a request that the store of uploaded inputs no longer holds.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct StoredInputs {
  /// Paths of the request that name a stored input whose file is gone, because no run references its contents any
  /// more, in the order of the request. Paths outside the store, such as files on the computer, are not listed.
  pub missing: Vec<PathBuf>,
}

/// Key of a file in the input store: the SHA-256 of its contents and its file name.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, JsonSchema, Serialize, Deserialize)]
pub struct StoredInput {
  /// SHA-256 of the file contents, as lowercase hexadecimal.
  pub sha256: String,
  /// File name of the upload.
  pub name: String,
}

pub struct UploadWriter {
  file: NamedTempFile,
  hasher: Sha256,
  name: String,
  size: usize,
  allowed: usize,
  limit: usize,
}

impl UploadWriter {
  pub fn write(&mut self, chunk: &[u8]) -> Result<(), Report> {
    self.size += chunk.len();
    if self.size > self.allowed {
      return Err(upload_too_large(self.limit, &self.name));
    }
    self.hasher.update(chunk);
    self
      .file
      .as_file_mut()
      .write_all(chunk)
      .wrap_err("When storing the upload")
  }

  pub fn finish(self) -> Result<StoredUpload, Report> {
    let Self {
      file,
      hasher,
      name,
      size,
      ..
    } = self;
    Ok(StoredUpload {
      file,
      key: StoredInput {
        sha256: hex_digest(&hasher.finalize())?,
        name,
      },
      size,
    })
  }
}

pub struct StoredUpload {
  file: NamedTempFile,
  key: StoredInput,
  size: usize,
}

impl StoredUpload {
  pub fn size(&self) -> usize {
    self.size
  }

  pub fn sha256(&self) -> &str {
    &self.key.sha256
  }
}

pub fn validate_upload_name(name: &str) -> Result<(), Report> {
  let mut components = Path::new(name).components();
  let single = matches!(components.next(), Some(Component::Normal(_))) && components.next().is_none();
  if !single || name.starts_with('.') || name.contains(['/', '\\']) {
    return Err(invalid(format!(
      "upload name `{name}` must be a plain file name without directories"
    )));
  }
  Ok(())
}

pub fn upload_too_large(limit: usize, name: &str) -> Report {
  Report::new(UploadTooLarge {
    limit: format_size(limit),
    name: name.to_owned(),
  })
}

pub fn format_size(bytes: usize) -> String {
  const UNITS: [&str; 5] = ["bytes", "KiB", "MiB", "GiB", "TiB"];
  let (value, unit) = UNITS
    .iter()
    .skip(1)
    .fold((bytes.to_f64_exact(), UNITS[0]), |(value, unit), next| {
      if value >= 1024.0 {
        (value / 1024.0, *next)
      } else {
        (value, unit)
      }
    });
  if unit == UNITS[0] {
    format!("{bytes} bytes")
  } else {
    format!("{} {unit}", float_to_significant_digits(value, 3))
  }
}

pub fn set_sha256<'a>(files: impl IntoIterator<Item = (&'a Path, &'a str)>) -> Result<String, Report> {
  let listing = files
    .into_iter()
    .map(|(path, sha256)| format!("{sha256}  {}\n", relative_name(path)))
    .sorted()
    .join("");
  Ok(hex_digest(&Sha256::digest(listing.as_bytes()))?)
}

fn is_sha256_hex(text: &str) -> bool {
  text.len() == SHA256_HEX_LEN && text.chars().all(|c| matches!(c, '0'..='9' | 'a'..='f'))
}

fn read_dir_if_exists(dir: &Path) -> Result<Vec<PathBuf>, Report> {
  let entries = match fs::read_dir(dir) {
    Ok(entries) => entries,
    Err(err) if err.kind() == io::ErrorKind::NotFound => return Ok(vec![]),
    Err(err) => return Err(Report::new(err).wrap_err(format!("When listing '{}'", dir.display()))),
  };
  entries
    .map(|entry| {
      entry
        .map(|entry| entry.path())
        .wrap_err_with(|| format!("When listing '{}'", dir.display()))
    })
    .collect()
}
