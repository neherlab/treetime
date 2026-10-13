use crate::runs::app_events::{AppChange, DATASETS_PATH, EXAMPLES_DOWNLOAD_PATH, StalePath};
use crate::runs::errors::conflict;
use crate::runs::manager::RunManager;
use app_commands::examples_download::{download_examples, ensure_empty_target};
use deser::Serialize;
use eyre::{Report, WrapErr};
use parking_lot::Mutex;
use schemars::JsonSchema;
use std::path::PathBuf;
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};
use treetime_schema::skip_serializing_optionals;
use treetime_utils::error::report_to_string;

const PROGRESS_INTERVAL: Duration = Duration::from_millis(250);

const DOWNLOAD_THREAD_NAME: &str = "treetime-examples";

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

  #[cfg_attr(
    dylint_lib = "custom",
    expect(
      spawn_handle_dropped,
      reason = "the thread publishes the download result as the download state that the app polls, and `download_examples` returns its failures as errors"
    )
  )]
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

/// The download of the example datasets, with the app event that reported it last.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct ExamplesDownloadStatus {
  pub download: ExamplesDownload,
  /// Sequence number of the last app event about the download. A client keeps whichever of this answer and the
  /// events it received has the higher number.
  pub seq: Option<u64>,
}

/// Stage of the download of the example datasets.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(tag = "state", rename_all = "kebab-case")]
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
