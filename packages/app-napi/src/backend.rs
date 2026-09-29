use crate::port::{PortHeader, PortReply, PortRequest};
use app_commands::bridge::error::ErrorResponse;
use app_commands::bridge::service::{AppService, Unconfined};
use app_commands::job::JobId;
use app_commands::runs::files::write_run_zip;
use app_commands::runs::manager::RunManager;
use app_server::routes::api_router;
use app_server::state::{DEFAULT_MAX_UPLOAD_SIZE, ServerConfig};
use axum::Router;
use eyre::{Report, WrapErr};
use napi::bindgen_prelude::Uint8Array;
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use tempfile::NamedTempFile;
use tokio::runtime::{Builder, Runtime};
use tokio::task::AbortHandle;
use tokio_stream::StreamExt as _;
use tokio_util::sync::CancellationToken;
use tower::ServiceExt as _;
use treetime_utils::env::env_var_optional;

const DATA_DIR_ENV: &str = "DATA_DIR";

const DEFAULT_DATA_DIR: &str = "data";

const RUNTIME_THREAD_NAME: &str = "treetime-desktop";

pub struct DesktopService {
  app: Arc<AppService>,
  router: Router,
  runtime: Runtime,
}

impl DesktopService {
  pub fn open(runs_dir: &Path) -> Result<Self, Report> {
    let data_dir = PathBuf::from(env_var_optional(DATA_DIR_ENV)?.unwrap_or_else(|| DEFAULT_DATA_DIR.to_owned()));
    let app = Arc::new(AppService::new(
      RunManager::open(runs_dir)?,
      data_dir.clone(),
      Arc::new(Unconfined),
    ));
    let (router, _) = api_router(
      Arc::clone(&app),
      ServerConfig {
        data_dir,
        runs_dir: runs_dir.to_path_buf(),
        max_upload_size: DEFAULT_MAX_UPLOAD_SIZE,
        shutdown: CancellationToken::new(),
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

fn write_atomically(destination: &Path, write: impl FnOnce(&mut File) -> Result<(), Report>) -> Result<(), Report> {
  let dir = destination
    .parent()
    .filter(|dir| !dir.as_os_str().is_empty())
    .unwrap_or_else(|| Path::new("."));
  let mut file = NamedTempFile::new_in(dir).wrap_err_with(|| format!("When creating a file in '{}'", dir.display()))?;
  write(file.as_file_mut())?;
  file
    .persist(destination)
    .wrap_err_with(|| format!("When saving '{}'", destination.display()))?;
  Ok(())
}
