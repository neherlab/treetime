use crate::runs::app_events::{AppChange, DATASETS_PATH, EXAMPLES_DOWNLOAD_PATH, StalePath};
use crate::runs::errors::conflict;
use crate::runs::manager::RunManager;
use crate::version::{BUILD_MODE, LONG_VERSION, NIGHTLY_BUILD_MODE};
use eyre::{Report, WrapErr};
use parking_lot::Mutex;
use reqwest::blocking::Client;
use rustls::crypto::CryptoProvider;
use rustls::crypto::ring::default_provider;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_with::skip_serializing_none;
use std::fs::{self, File};
use std::io::{ErrorKind, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};
use tempfile::Builder;
use treetime_utils::error::report_to_string;
use treetime_utils::{make_error, make_report};
use zip::ZipArchive;

pub const EXAMPLES_ASSET: &str = "examples.zip";

const RELEASES_URL: &str = "https://github.com/neherlab/treetime-nightly/releases";

const TEMPORARY_PREFIX: &str = ".treetime-examples-";

const FINDER_METADATA_FILE: &str = ".DS_Store";

const CHUNK_SIZE: usize = 1 << 16;

const PROGRESS_INTERVAL: Duration = Duration::from_millis(250);

const DOWNLOAD_THREAD_NAME: &str = "treetime-examples";

const CONNECT_TIMEOUT: Duration = Duration::from_secs(30);

pub fn examples_url() -> String {
  if BUILD_MODE == NIGHTLY_BUILD_MODE {
    format!(
      "{RELEASES_URL}/download/{}/{EXAMPLES_ASSET}",
      LONG_VERSION.replace('+', "%2B")
    )
  } else {
    format!("{RELEASES_URL}/latest/download/{EXAMPLES_ASSET}")
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DownloadProgress {
  pub received: u64,
  pub total: Option<u64>,
}

pub fn download_examples(url: &str, target: &Path, mut progress: impl FnMut(DownloadProgress)) -> Result<(), Report> {
  ensure_empty_target(target)?;
  let parent = target
    .parent()
    .filter(|parent| !parent.as_os_str().is_empty())
    .unwrap_or_else(|| Path::new("."));
  fs::create_dir_all(parent).wrap_err_with(|| format!("When creating the folder '{}'", parent.display()))?;
  remove_stale_temporaries(parent)?;

  let mut archive = Builder::new()
    .prefix(TEMPORARY_PREFIX)
    .suffix(".zip")
    .tempfile_in(parent)
    .wrap_err_with(|| format!("When creating a temporary file in '{}'", parent.display()))?;
  fetch(url, archive.as_file_mut(), &mut progress).wrap_err_with(|| format!("When downloading '{url}'"))?;
  archive.as_file_mut().seek(SeekFrom::Start(0))?;

  let unpacked = Builder::new()
    .prefix(TEMPORARY_PREFIX)
    .tempdir_in(parent)
    .wrap_err_with(|| format!("When creating a temporary folder in '{}'", parent.display()))?;
  ZipArchive::new(archive.as_file_mut())
    .and_then(|mut zip| zip.extract(unpacked.path()))
    .wrap_err_with(|| format!("When unpacking the example datasets from '{url}'"))?;

  remove_empty_target(target)?;
  fs::rename(unpacked.path(), target).wrap_err_with(|| {
    format!(
      "When moving the example datasets from '{}' to '{}'",
      unpacked.path().display(),
      target.display()
    )
  })?;
  let _moved = unpacked.keep();
  Ok(())
}

#[derive(Debug)]
pub struct ExampleDownloads {
  url: String,
  target: PathBuf,
  state: Mutex<ExamplesDownloadStatus>,
}

impl ExampleDownloads {
  pub fn new(url: impl Into<String>, target: impl Into<PathBuf>) -> Self {
    Self {
      url: url.into(),
      target: target.into(),
      state: Mutex::new(ExamplesDownloadStatus {
        download: ExamplesDownload::Idle,
        seq: None,
      }),
    }
  }

  pub fn status(&self) -> ExamplesDownloadStatus {
    self.state.lock().clone()
  }

  pub fn start(self: &Arc<Self>, runs: &Arc<RunManager>) -> Result<ExamplesDownloadStatus, Report> {
    let mut state = self.state.lock();
    if matches!(state.download, ExamplesDownload::Running { .. }) {
      return Err(conflict("the example datasets are being downloaded already"));
    }
    ensure_empty_target(&self.target).map_err(|report| conflict(report_to_string(&report)))?;
    publish(
      &mut state,
      runs,
      ExamplesDownload::Running {
        received: 0,
        total: None,
      },
    );
    let status = state.clone();
    drop(state);
    let downloads = Arc::clone(self);
    let runs = Arc::clone(runs);
    thread::Builder::new()
      .name(DOWNLOAD_THREAD_NAME.to_owned())
      .spawn(move || downloads.run(&runs))
      .wrap_err("When starting the download of the example datasets")?;
    Ok(status)
  }

  fn run(&self, runs: &Arc<RunManager>) {
    let mut reported = Instant::now();
    let mut received = 0;
    let result = download_examples(&self.url, &self.target, |progress| {
      received = progress.received;
      if reported.elapsed() >= PROGRESS_INTERVAL {
        reported = Instant::now();
        publish(
          &mut self.state.lock(),
          runs,
          ExamplesDownload::Running {
            received: progress.received,
            total: progress.total,
          },
        );
      }
    });
    let download = match result {
      Ok(()) => ExamplesDownload::Done { received },
      Err(report) => ExamplesDownload::Failed {
        message: report_to_string(&report),
      },
    };
    publish(&mut self.state.lock(), runs, download);
  }
}

/// Stage of the download of the example datasets.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
#[deser(skip_serializing_optionals)]
#[serde(tag = "state", rename_all = "kebab-case")]
#[deser(tag = "state", rename_all = "kebab-case")]
pub enum ExamplesDownload {
  /// No download started since the back end started.
  Idle,
  /// The archive is being downloaded or unpacked.
  Running {
    /// Bytes received so far.
    received: u64,
    /// Size of the archive in bytes, when the server sends it.
    total: Option<u64>,
  },
  /// The example datasets are in the examples folder.
  Done {
    /// Size of the archive in bytes.
    received: u64,
  },
  /// The download failed; the examples folder is unchanged.
  Failed {
    /// What failed.
    message: String,
  },
}

/// The download of the example datasets, with the app event that reported it last.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
#[deser(skip_serializing_optionals)]
pub struct ExamplesDownloadStatus {
  pub download: ExamplesDownload,
  /// Sequence number of the last app event about the download. A client keeps whichever of this answer and the
  /// events it received has the higher number.
  pub seq: Option<usize>,
}

fn publish(state: &mut ExamplesDownloadStatus, runs: &RunManager, download: ExamplesDownload) {
  let mut stale = vec![StalePath::exact(EXAMPLES_DOWNLOAD_PATH)];
  if matches!(download, ExamplesDownload::Done { .. }) {
    stale.push(StalePath::exact(DATASETS_PATH));
  }
  let event = runs.app_events().append(
    AppChange::ExamplesDownload {
      download: download.clone(),
    },
    stale,
  );
  *state = ExamplesDownloadStatus {
    download,
    seq: Some(event.seq),
  };
}

fn ensure_empty_target(target: &Path) -> Result<(), Report> {
  let entries = match fs::read_dir(target) {
    Ok(entries) => entries,
    Err(error) if error.kind() == ErrorKind::NotFound => return Ok(()),
    Err(error) => {
      return Err(Report::new(error)).wrap_err_with(|| format!("When reading the folder '{}'", target.display()));
    },
  };
  for entry in entries {
    let entry = entry.wrap_err_with(|| format!("When reading the folder '{}'", target.display()))?;
    if entry.file_name() != FINDER_METADATA_FILE {
      return make_error!(
        "the folder '{}' is not empty; the example datasets go into an empty or missing folder",
        target.display()
      );
    }
  }
  Ok(())
}

fn remove_stale_temporaries(parent: &Path) -> Result<(), Report> {
  let entries = fs::read_dir(parent).wrap_err_with(|| format!("When reading the folder '{}'", parent.display()))?;
  for entry in entries {
    let entry = entry.wrap_err_with(|| format!("When reading the folder '{}'", parent.display()))?;
    if !entry.file_name().to_string_lossy().starts_with(TEMPORARY_PREFIX) {
      continue;
    }
    let path = entry.path();
    let removed = if entry.file_type()?.is_dir() {
      fs::remove_dir_all(&path)
    } else {
      fs::remove_file(&path)
    };
    removed.wrap_err_with(|| format!("When removing the stale temporary '{}'", path.display()))?;
  }
  Ok(())
}

fn fetch(url: &str, file: &mut File, progress: &mut impl FnMut(DownloadProgress)) -> Result<(), Report> {
  install_crypto_provider()?;
  let client = Client::builder()
    .user_agent(concat!("treetime/", env!("CARGO_PKG_VERSION")))
    .connect_timeout(CONNECT_TIMEOUT)
    .timeout(None)
    .build()
    .wrap_err("When creating the HTTP client")?;
  let mut response = client.get(url).send()?.error_for_status()?;
  let total = response.content_length();
  let mut received = 0;
  let mut buffer = vec![0; CHUNK_SIZE];
  progress(DownloadProgress { received, total });
  loop {
    let read = response.read(&mut buffer)?;
    if read == 0 {
      break;
    }
    file.write_all(&buffer[..read])?;
    received += u64::try_from(read)?;
    progress(DownloadProgress { received, total });
  }
  file.flush()?;
  if let Some(total) = total.filter(|total| *total != received) {
    return make_error!("the download ended after {received} of {total} bytes");
  }
  Ok(())
}

fn install_crypto_provider() -> Result<(), Report> {
  if CryptoProvider::get_default().is_some() || default_provider().install_default().is_ok() {
    return Ok(());
  }
  CryptoProvider::get_default()
    .map(|_| ())
    .ok_or_else(|| make_report!("the TLS crypto provider could not be installed"))
}

fn remove_empty_target(target: &Path) -> Result<(), Report> {
  match fs::remove_file(target.join(FINDER_METADATA_FILE)) {
    Ok(()) => {},
    Err(error) if error.kind() == ErrorKind::NotFound => {},
    Err(error) => {
      return Err(Report::new(error)).wrap_err_with(|| format!("When cleaning the folder '{}'", target.display()));
    },
  }
  match fs::remove_dir(target) {
    Ok(()) => Ok(()),
    Err(error) if error.kind() == ErrorKind::NotFound => Ok(()),
    Err(error) => Err(Report::new(error)).wrap_err_with(|| {
      format!(
        "When replacing the folder '{}', which is no longer empty",
        target.display()
      )
    }),
  }
}
