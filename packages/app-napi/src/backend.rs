use crate::port::{PortHeader, PortReply, PortRequest};
use app_commands::app_paths::{AppFolderEnv, AppFolderPath, AppPaths};
use app_commands::app_settings::store::AppSettingsStore;
use app_commands::app_settings::workspace::{OpenedRuns, open_runs_folder};
use app_commands::atomic_write::write_atomically;
use app_commands::bridge::error::ErrorResponse;
use app_commands::bridge::service::{AppService, Unconfined};
use app_commands::job::JobId;
use app_commands::runs::files::write_run_zip;
use app_server::routes::api_router;
use app_server::state::{DEFAULT_MAX_UPLOAD_SIZE, LocalSettings, ServerConfig};
use axum::Router;
use eyre::{Report, WrapErr};
use napi::bindgen_prelude::Uint8Array;
use std::fs::{self, File};
use std::io;
use std::path::Path;
use std::sync::Arc;
use tokio::runtime::{Builder, Runtime};
use tokio::task::AbortHandle;
use tokio_stream::StreamExt as _;
use tokio_util::sync::CancellationToken;
use tower::ServiceExt as _;

const RUNTIME_THREAD_NAME: &str = "treetime-desktop";

pub struct DesktopService {
  app: Arc<AppService>,
  router: Router,
  runtime: Runtime,
}

impl DesktopService {
  pub fn open(root: &Path, env: &AppFolderEnv) -> Result<Self, Report> {
    let store = Arc::new(AppSettingsStore::open(root)?);
    let settings = store.read()?;
    let mut paths = AppPaths::resolve(root, env, &settings.paths);
    let named_in = settings.paths.runs.is_some().then(|| store.path());
    let OpenedRuns {
      runs,
      error: runs_error,
    } = open_runs_folder(&mut paths, named_in)?;
    fs::create_dir_all(&paths.examples.path).wrap_err_with(|| {
      format!(
        "When creating the examples folder '{}'{}",
        paths.examples.path.display(),
        folder_source(&paths.examples, settings.paths.examples.is_some(), store.path())
      )
    })?;
    let examples = paths.examples.path.clone();
    let app = Arc::new(AppService::new(runs, examples.clone(), Arc::new(Unconfined)));
    let (router, _) = api_router(
      Arc::clone(&app),
      ServerConfig {
        examples_dir: examples,
        runs_dir: paths.runs.path.clone(),
        max_upload_size: DEFAULT_MAX_UPLOAD_SIZE,
        shutdown: CancellationToken::new(),
        settings: Some(LocalSettings {
          store,
          paths,
          runs_error,
        }),
      },
    )?;
    let runtime = Builder::new_multi_thread()
      .enable_all()
      .thread_name(RUNTIME_THREAD_NAME)
      .build()
      .wrap_err("When starting the runtime of the desktop back end")?;
    Ok(Self { app, router, runtime })
  }

  pub fn fetch(&self, request: PortRequest, send: impl Fn(PortReply) -> bool + Send + Sync + 'static) -> AbortHandle {
    let seq = request.seq;
    let send = Arc::new(send);
    let exchange = self
      .runtime
      .spawn(exchange(self.router.clone(), request, Arc::clone(&send)));
    let abort = exchange.abort_handle();
    self.runtime.spawn(async move {
      if let Err(err) = exchange.await
        && err.is_panic()
      {
        send(PortReply::error(
          seq,
          ErrorResponse::from_panic(err.into_panic().as_ref()),
        ));
      }
    });
    abort
  }

  pub fn save_run_file(&self, id: &JobId, relative: &str, destination: &Path) -> Result<(), Report> {
    let source = self.app.runs().file_path(id, relative)?;
    write_atomically(destination, |file| {
      let mut reader = File::open(&source).wrap_err_with(|| format!("When opening '{}'", source.display()))?;
      io::copy(&mut reader, file)?;
      Ok(())
    })
  }

  pub fn save_run_archive(&self, id: &JobId, destination: &Path) -> Result<(), Report> {
    let runs = self.app.runs();
    runs.get(id)?;
    let out_dir = runs.store().out_dir(id);
    write_atomically(destination, |file| write_run_zip(&out_dir, id.as_str(), file))
  }
}

async fn exchange(router: Router, request: PortRequest, send: Arc<impl Fn(PortReply) -> bool + Send + Sync>) {
  let seq = request.seq;
  let request = match request.into_http() {
    Ok(request) => request,
    Err(report) => {
      send(PortReply::error(seq, ErrorResponse::from_report(&report)));
      return;
    },
  };
  let response = match router.oneshot(request).await {
    Ok(response) => response,
    Err(infallible) => match infallible {},
  };
  let (head, body) = response.into_parts();
  let head = PortReply::Head {
    seq,
    status: head.status.as_u16(),
    headers: PortHeader::list(&head.headers),
  };
  if !send(head) {
    return;
  }
  let mut frames = body.into_data_stream();
  while let Some(frame) = frames.next().await {
    let reply = match frame {
      Ok(bytes) => PortReply::Chunk {
        seq,
        data: Uint8Array::new(Vec::from(bytes)),
      },
      Err(err) => {
        let report = Report::new(err).wrap_err("When reading the response body");
        send(PortReply::error(seq, ErrorResponse::from_report(&report)));
        return;
      },
    };
    if !send(reply) {
      return;
    }
  }
  send(PortReply::End { seq });
}

fn folder_source(folder: &AppFolderPath, in_settings: bool, settings_file: &Path) -> String {
  match folder.fixed_by {
    Some(variable) => format!(", set by {variable}"),
    None if in_settings => format!(", named in '{}'", settings_file.display()),
    None => String::new(),
  }
}
